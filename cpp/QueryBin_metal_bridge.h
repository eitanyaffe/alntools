#ifndef QUERYBIN_METAL_BRIDGE_H
#define QUERYBIN_METAL_BRIDGE_H

#include <stdint.h>

// C-compatible structs for GPU input
typedef struct {
  uint32_t contig_index;
  uint32_t contig_start;
  uint32_t contig_end;
  uint32_t num_mutations;
} MetalSimpleAlignment;

typedef struct {
  uint32_t start;
  uint32_t end;
} MetalSimpleInterval;

#ifdef __cplusplus
extern "C" {
#endif

// Forward declaration
typedef struct MetalBinImplementation MetalBinImplementation;

// C interface for Metal implementation
MetalBinImplementation* metal_create_implementation(void);
void metal_destroy_implementation(MetalBinImplementation* impl);

int metal_aggregate_data(MetalBinImplementation* impl,
                        const MetalSimpleAlignment* alignments,
                        int num_alignments,
                        const MetalSimpleInterval* intervals,
                        int num_intervals,
                        const uint32_t* contig_indices,
                        uint32_t binsize,
                        uint32_t max_contigs,
                        uint32_t max_bins_per_contig);

// Result extraction functions
int metal_get_result_count(MetalBinImplementation* impl);
void metal_get_result(MetalBinImplementation* impl, int index,
                     uint32_t* contig_index, uint32_t* bin_start,
                     uint32_t* sequenced_bp, uint32_t* mutation_count,
                     uint32_t* dist_none, uint32_t* dist_5, uint32_t* dist_4,
                     uint32_t* dist_3, uint32_t* dist_2, uint32_t* dist_1_plus);

#ifdef __cplusplus
}
#endif

#endif // QUERYBIN_METAL_BRIDGE_H
