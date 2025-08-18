#ifndef QUERYBIN_GPU_H
#define QUERYBIN_GPU_H

#include "QueryBin.h"
#include <memory>

#ifdef METAL_SUPPORT
    // Forward declaration for Metal implementation
    struct MetalBinImplementation;
#endif

// GPU-accelerated version of QueryBin - automatically uses best available GPU
class QueryBinGPU {
private:
    std::unique_ptr<QueryBin> cpu_fallback;
    
    bool use_gpu;
    
#ifdef METAL_SUPPORT
    MetalBinImplementation* metal_impl;
#else
    void* metal_impl;  // placeholder when not compiled with Metal support
#endif

    const std::vector<Interval>& intervals;
    const AlignmentStore& store;
    int binsize;
    double seg_threshold;
    double non_ref_threshold;
    ClipMode clip_mode;
    int clip_margin;
    double min_mutations_percent;
    double max_mutations_percent;
    
    std::map<std::pair<uint32_t, uint32_t>, BinData> bin_results;
    std::vector<BinOutputRow> output_rows;

public:
    QueryBinGPU(
        const std::vector<Interval>& intervals,
        const AlignmentStore& store,
        int binsize,
        double seg_threshold = 0.2,
        double non_ref_threshold = 0.9,
        int num_threads = 0,  // ignored for GPU version
        ClipMode clip_mode = ClipMode::ALL,
        int clip_margin = 10,
        double min_mutations_percent = 0.0,
        double max_mutations_percent = 10.0);
    
    ~QueryBinGPU();
    
    // main interface - same as QueryBin
    void execute();
    void write_to_csv(const std::string& ofn_prefix);
    const std::vector<BinOutputRow>& get_output_rows() const { return output_rows; }
    
private:
    void aggregate_data();
    void generate_output_rows();
};

#endif // QUERYBIN_GPU_H
