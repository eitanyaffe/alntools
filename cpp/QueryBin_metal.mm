#include "QueryBin_metal_bridge.h"
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>
#import <Foundation/Foundation.h>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <vector>

// Result storage
struct MetalBinResult {
    uint32_t contig_index;
    uint32_t bin_start;
    uint32_t sequenced_bp;
    uint32_t mutation_count;
    uint32_t dist_none;
    uint32_t dist_5;
    uint32_t dist_4;
    uint32_t dist_3;
    uint32_t dist_2;
    uint32_t dist_1_plus;
};

// Metal implementation class
@interface MetalBinImplementationObjC : NSObject {
    id<MTLDevice> device;
    id<MTLCommandQueue> commandQueue;
    id<MTLComputePipelineState> pipelineState;
    id<MTLLibrary> library;
    std::vector<MetalBinResult> results;
}

- (instancetype)init;
- (BOOL)setupComputeShader;
- (BOOL)processAlignments:(const MetalSimpleAlignment*)alignments
                    count:(int)num_alignments
                intervals:(const MetalSimpleInterval*)intervals
            intervalCount:(int)num_intervals
            contigIndices:(const uint32_t*)contig_indices
                  binsize:(uint32_t)binsize
               maxContigs:(uint32_t)max_contigs
        maxBinsPerContig:(uint32_t)max_bins_per_contig;
- (int)getResultCount;
- (MetalBinResult)getResultAtIndex:(int)index;

@end

@implementation MetalBinImplementationObjC

- (instancetype)init {
    self = [super init];
    if (self) {
        device = MTLCreateSystemDefaultDevice();
        if (!device) {
            NSLog(@"Metal device not available");
            return nil;
        }
        
        commandQueue = [device newCommandQueue];
        
        NSLog(@"Initialized Metal compute on Apple Silicon GPU: %@", [device name]);
        
        if (![self setupComputeShader]) {
            return nil;
        }
    }
    return self;
}

- (BOOL)setupComputeShader {
    // Simplified Metal shader for now - just does basic bin counting
    NSString* shaderSource = @R"(
#include <metal_stdlib>
using namespace metal;

struct SimpleAlignment { uint32_t contig_index, contig_start, contig_end, num_mutations; };
struct SimpleInterval { uint32_t start, end; };

kernel void process_alignments_simple(
    device const SimpleAlignment* alignments [[buffer(0)]],
    device const SimpleInterval* intervals [[buffer(1)]],
    device const uint32_t* contig_indices [[buffer(2)]],
    device atomic_uint* sequenced_bp [[buffer(3)]],
    device atomic_uint* mutation_count [[buffer(4)]],
    device atomic_uint* dist_none [[buffer(5)]],
    device atomic_uint* dist_5 [[buffer(6)]],
    device atomic_uint* dist_4 [[buffer(7)]],
    device atomic_uint* dist_3 [[buffer(8)]],
    device atomic_uint* dist_2 [[buffer(9)]],
    device atomic_uint* dist_1_plus [[buffer(10)]],
    constant uint32_t& binsize [[buffer(11)]],
    constant uint32_t& num_intervals [[buffer(12)]],
    constant uint32_t& max_bins_per_contig [[buffer(13)]],
    constant uint32_t& max_contigs [[buffer(14)]],
    constant uint32_t& num_alignments [[buffer(15)]],
    uint tid [[thread_position_in_grid]]
) {
    if (tid >= num_alignments) return;
    
    const SimpleAlignment aln = alignments[tid];
    
    // process each interval
    for (uint32_t i = 0; i < num_intervals; i++) {
        const SimpleInterval interval = intervals[i];
        
        // skip if different contig (simplified check)
        if (aln.contig_index != contig_indices[i]) continue;
        
        // skip if no overlap
        if (aln.contig_end <= interval.start || aln.contig_start >= interval.end) continue;
        
        uint32_t adjusted_start = (interval.start / binsize) * binsize;
        uint32_t last_bin_start = ((interval.end - 1) / binsize) * binsize;
        
        // process bins
        for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += binsize) {
            uint32_t b_end = b_start + binsize;
            
            // calculate overlap
            uint32_t effective_start = max(max(aln.contig_start, b_start), interval.start);
            uint32_t effective_end = min(min(aln.contig_end, b_end), interval.end);
            
            if (effective_end > effective_start) {
                int overlap_length = effective_end - effective_start;
                
                // calculate bin index
                uint32_t bin_id_in_contig = (b_start / binsize);
                if (aln.contig_index < max_contigs && bin_id_in_contig < max_bins_per_contig) {
                    uint32_t bin_index = aln.contig_index * max_bins_per_contig + bin_id_in_contig;
                    // atomic updates
                    atomic_fetch_add_explicit(&sequenced_bp[bin_index], overlap_length, memory_order_relaxed);
                    atomic_fetch_add_explicit(&mutation_count[bin_index], aln.num_mutations, memory_order_relaxed);
                    
                    // distance categories by mutations per bp
                    float aln_len = (float)max((uint32_t)1, aln.contig_end - aln.contig_start);
                    float mpbp = (float)aln.num_mutations / aln_len;
                    if (aln.num_mutations == 0) {
                        atomic_fetch_add_explicit(&dist_none[bin_index], 1, memory_order_relaxed);
                    } else if (mpbp < 1e-4f) {
                        atomic_fetch_add_explicit(&dist_5[bin_index], 1, memory_order_relaxed);
                    } else if (mpbp < 1e-3f) {
                        atomic_fetch_add_explicit(&dist_4[bin_index], 1, memory_order_relaxed);
                    } else if (mpbp < 1e-2f) {
                        atomic_fetch_add_explicit(&dist_3[bin_index], 1, memory_order_relaxed);
                    } else if (mpbp < 1e-1f) {
                        atomic_fetch_add_explicit(&dist_2[bin_index], 1, memory_order_relaxed);
                    } else {
                        atomic_fetch_add_explicit(&dist_1_plus[bin_index], 1, memory_order_relaxed);
                    }
                }
            }
        }
    }
}
)";
    
    NSError* error = nil;
    library = [device newLibraryWithSource:shaderSource options:nil error:&error];
    
    if (!library) {
        NSLog(@"Error creating Metal library: %@", [error localizedDescription]);
        return NO;
    }
    
    id<MTLFunction> kernelFunction = [library newFunctionWithName:@"process_alignments_simple"];
    if (!kernelFunction) {
        NSLog(@"Could not find kernel function");
        return NO;
    }
    
    pipelineState = [device newComputePipelineStateWithFunction:kernelFunction error:&error];
    
    if (!pipelineState) {
        NSLog(@"Error creating compute pipeline: %@", [error localizedDescription]);
        return NO;
    }
    
    NSLog(@"Metal compute shader compiled successfully");
    return YES;
}

- (BOOL)processAlignments:(const MetalSimpleAlignment*)alignments
                    count:(int)num_alignments
                intervals:(const MetalSimpleInterval*)intervals
            intervalCount:(int)num_intervals
            contigIndices:(const uint32_t*)contig_indices
                  binsize:(uint32_t)binsize
               maxContigs:(uint32_t)max_contigs
        maxBinsPerContig:(uint32_t)max_bins_per_contig {
    
    // Use full requested sizes
    int safe_alignment_count = num_alignments;
    int safe_interval_count = num_intervals;
    uint32_t safe_bins = max_bins_per_contig;
    size_t total_bins = (size_t)max_contigs * (size_t)safe_bins;
    
    NSLog(@"Processing %d alignments, %d intervals, %zu total bins", 
          safe_alignment_count, safe_interval_count, total_bins);
    
    // Convert alignments to simplified format
    struct SimpleAlignment {
        uint32_t contig_index;
        uint32_t contig_start;
        uint32_t contig_end;
        uint32_t num_mutations;
    };
    
    struct SimpleInterval {
        uint32_t start;
        uint32_t end;
    };
    
    // Copy directly; caller provides simplified structs
    std::vector<SimpleAlignment> simple_alignments;
    simple_alignments.reserve(safe_alignment_count);
    for (int i = 0; i < safe_alignment_count; i++) {
        SimpleAlignment s = {alignments[i].contig_index, alignments[i].contig_start, alignments[i].contig_end, alignments[i].num_mutations};
        simple_alignments.push_back(s);
    }
    std::vector<SimpleInterval> simple_intervals;
    simple_intervals.reserve(safe_interval_count);
    for (int i = 0; i < safe_interval_count; i++) simple_intervals.push_back({intervals[i].start, intervals[i].end});
    
    // Create Metal buffers
    id<MTLBuffer> alignmentBuffer = [device newBufferWithBytes:simple_alignments.data()
                                                        length:safe_alignment_count * sizeof(SimpleAlignment)
                                                       options:MTLResourceStorageModeShared];
    
    id<MTLBuffer> intervalBuffer = [device newBufferWithBytes:simple_intervals.data()
                                                       length:safe_interval_count * sizeof(SimpleInterval)
                                                      options:MTLResourceStorageModeShared];
    
    id<MTLBuffer> contigIndicesBuffer = [device newBufferWithBytes:contig_indices
                                                            length:safe_interval_count * sizeof(uint32_t)
                                                           options:MTLResourceStorageModeShared];
    
    // Create result buffers
    id<MTLBuffer> sequencedBpBuffer = [device newBufferWithLength:total_bins * sizeof(uint32_t)
                                                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> mutationCountBuffer = [device newBufferWithLength:total_bins * sizeof(uint32_t)
                                                            options:MTLResourceStorageModeShared];
    id<MTLBuffer> distNoneBuffer = [device newBufferWithLength:total_bins * sizeof(uint32_t)
                                                       options:MTLResourceStorageModeShared];
    id<MTLBuffer> dist5Buffer = [device newBufferWithLength:total_bins * sizeof(uint32_t)
                                                     options:MTLResourceStorageModeShared];
    id<MTLBuffer> dist4Buffer = [device newBufferWithLength:total_bins * sizeof(uint32_t)
                                                     options:MTLResourceStorageModeShared];
    id<MTLBuffer> dist3Buffer = [device newBufferWithLength:total_bins * sizeof(uint32_t)
                                                     options:MTLResourceStorageModeShared];
    id<MTLBuffer> dist2Buffer = [device newBufferWithLength:total_bins * sizeof(uint32_t)
                                                     options:MTLResourceStorageModeShared];
    id<MTLBuffer> dist1PlusBuffer = [device newBufferWithLength:total_bins * sizeof(uint32_t)
                                                          options:MTLResourceStorageModeShared];
    
    // Zero out result buffers
    memset([sequencedBpBuffer contents], 0, total_bins * sizeof(uint32_t));
    memset([mutationCountBuffer contents], 0, total_bins * sizeof(uint32_t));
    memset([distNoneBuffer contents], 0, total_bins * sizeof(uint32_t));
    memset([dist5Buffer contents], 0, total_bins * sizeof(uint32_t));
    memset([dist4Buffer contents], 0, total_bins * sizeof(uint32_t));
    memset([dist3Buffer contents], 0, total_bins * sizeof(uint32_t));
    memset([dist2Buffer contents], 0, total_bins * sizeof(uint32_t));
    memset([dist1PlusBuffer contents], 0, total_bins * sizeof(uint32_t));
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Create command buffer and encoder
    id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    
    [encoder setComputePipelineState:pipelineState];
    
    // Set buffers
    [encoder setBuffer:alignmentBuffer offset:0 atIndex:0];
    [encoder setBuffer:intervalBuffer offset:0 atIndex:1];
    [encoder setBuffer:contigIndicesBuffer offset:0 atIndex:2];
    [encoder setBuffer:sequencedBpBuffer offset:0 atIndex:3];
    [encoder setBuffer:mutationCountBuffer offset:0 atIndex:4];
    [encoder setBuffer:distNoneBuffer offset:0 atIndex:5];
    [encoder setBuffer:dist5Buffer offset:0 atIndex:6];
    [encoder setBuffer:dist4Buffer offset:0 atIndex:7];
    [encoder setBuffer:dist3Buffer offset:0 atIndex:8];
    [encoder setBuffer:dist2Buffer offset:0 atIndex:9];
    [encoder setBuffer:dist1PlusBuffer offset:0 atIndex:10];
    
    // Set constants
    [encoder setBytes:&binsize length:sizeof(uint32_t) atIndex:11];
    uint32_t num_intervals_u32 = safe_interval_count;
    [encoder setBytes:&num_intervals_u32 length:sizeof(uint32_t) atIndex:12];
    [encoder setBytes:&safe_bins length:sizeof(uint32_t) atIndex:13];
    [encoder setBytes:&max_contigs length:sizeof(uint32_t) atIndex:14];
    uint32_t num_alignments_u32 = safe_alignment_count;
    [encoder setBytes:&num_alignments_u32 length:sizeof(uint32_t) atIndex:15];
    
    // Dispatch exact number of threads for alignments
    MTLSize threadsPerGroup = MTLSizeMake(64, 1, 1);
    MTLSize threads = MTLSizeMake((NSUInteger)safe_alignment_count, 1, 1);
    
    NSLog(@"Dispatching %lu threads with %lu threads per group",
          (unsigned long)threads.width, (unsigned long)threadsPerGroup.width);
    
    [encoder dispatchThreads:threads threadsPerThreadgroup:threadsPerGroup];
    [encoder endEncoding];
    
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    NSLog(@"Metal GPU processing completed in %lld ms", duration.count());
    
    if (commandBuffer.status != MTLCommandBufferStatusCompleted) {
        NSLog(@"Metal command buffer execution failed");
        return NO;
    }
    
    // Extract results
    results.clear();
    uint32_t* sequenced_bp = (uint32_t*)[sequencedBpBuffer contents];
    uint32_t* mutation_count = (uint32_t*)[mutationCountBuffer contents];
    uint32_t* dist_none = (uint32_t*)[distNoneBuffer contents];
    uint32_t* dist_5 = (uint32_t*)[dist5Buffer contents];
    uint32_t* dist_4 = (uint32_t*)[dist4Buffer contents];
    uint32_t* dist_3 = (uint32_t*)[dist3Buffer contents];
    uint32_t* dist_2 = (uint32_t*)[dist2Buffer contents];
    uint32_t* dist_1_plus = (uint32_t*)[dist1PlusBuffer contents];
    
    for (uint32_t contig_idx = 0; contig_idx < max_contigs; contig_idx++) {
        for (uint32_t bin_idx = 0; bin_idx < safe_bins; bin_idx++) {
            uint32_t linear_idx = contig_idx * safe_bins + bin_idx;
            
            if (linear_idx < total_bins && sequenced_bp[linear_idx] > 0) {
                MetalBinResult result = {
                    contig_idx,
                    bin_idx * binsize,
                    sequenced_bp[linear_idx],
                    mutation_count[linear_idx],
                    dist_none[linear_idx],
                    dist_5[linear_idx],
                    dist_4[linear_idx],
                    dist_3[linear_idx],
                    dist_2[linear_idx],
                    dist_1_plus[linear_idx]
                };
                results.push_back(result);
            }
        }
    }
    
    NSLog(@"Extracted %zu non-empty bins from GPU results", results.size());
    return YES;
}

- (int)getResultCount {
    return static_cast<int>(results.size());
}

- (MetalBinResult)getResultAtIndex:(int)index {
    if (index >= 0 && index < static_cast<int>(results.size())) {
        return results[index];
    }
    MetalBinResult empty = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    return empty;
}

@end

// C interface implementation
struct MetalBinImplementation {
    MetalBinImplementationObjC* objc_impl;
};

MetalBinImplementation* metal_create_implementation(void) {
    MetalBinImplementation* impl = new MetalBinImplementation();
    impl->objc_impl = [[MetalBinImplementationObjC alloc] init];
    
    if (!impl->objc_impl) {
        delete impl;
        return nullptr;
    }
    
    return impl;
}

void metal_destroy_implementation(MetalBinImplementation* impl) {
    if (impl) {
        impl->objc_impl = nil;
        delete impl;
    }
}

int metal_aggregate_data(MetalBinImplementation* impl,
                        const MetalSimpleAlignment* alignments,
                        int num_alignments,
                        const MetalSimpleInterval* intervals,
                        int num_intervals,
                        const uint32_t* contig_indices,
                        uint32_t binsize,
                        uint32_t max_contigs,
                        uint32_t max_bins_per_contig) {
    
    if (!impl || !impl->objc_impl) return 0;
    
    BOOL success = [impl->objc_impl processAlignments:alignments
                                                count:num_alignments
                                            intervals:intervals
                                        intervalCount:num_intervals
                                        contigIndices:contig_indices
                                              binsize:binsize
                                           maxContigs:max_contigs
                                     maxBinsPerContig:max_bins_per_contig];
    
    return success ? 1 : 0;
}

int metal_get_result_count(MetalBinImplementation* impl) {
    if (!impl || !impl->objc_impl) return 0;
    return [impl->objc_impl getResultCount];
}

void metal_get_result(MetalBinImplementation* impl, int index,
                     uint32_t* contig_index, uint32_t* bin_start,
                     uint32_t* sequenced_bp, uint32_t* mutation_count,
                     uint32_t* dist_none, uint32_t* dist_5, uint32_t* dist_4,
                     uint32_t* dist_3, uint32_t* dist_2, uint32_t* dist_1_plus) {
    
    if (!impl || !impl->objc_impl) return;
    
    MetalBinResult result = [impl->objc_impl getResultAtIndex:index];
    
    if (contig_index) *contig_index = result.contig_index;
    if (bin_start) *bin_start = result.bin_start;
    if (sequenced_bp) *sequenced_bp = result.sequenced_bp;
    if (mutation_count) *mutation_count = result.mutation_count;
    if (dist_none) *dist_none = result.dist_none;
    if (dist_5) *dist_5 = result.dist_5;
    if (dist_4) *dist_4 = result.dist_4;
    if (dist_3) *dist_3 = result.dist_3;
    if (dist_2) *dist_2 = result.dist_2;
    if (dist_1_plus) *dist_1_plus = result.dist_1_plus;
}