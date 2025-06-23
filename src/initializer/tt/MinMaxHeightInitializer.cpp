#include "initializer/tt/MinMaxHeightInitializer.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(MinMaxHeightInitializer);

bool MinMaxHeightInitializer::operator()(const std::vector<vec3>& hits) {
    std::size_t min_idx = 0, max_idx = 0;

    for (std::size_t k = 1; k < hits.size(); ++k) {
        if (hits[k].z < hits[min_idx].z) min_idx = k;
        if (hits[max_idx].z < hits[k].z) max_idx = k;
    }

    if (min_idx == max_idx) {
        LogError << "Min and max hits are the same\n";
        return false;
    }

    m_params = std::vector<double>{
        hits[max_idx].x, 
        hits[max_idx].y, 
        hits[max_idx].z,
        hits[min_idx].x - hits[max_idx].x, 
        hits[min_idx].y - hits[max_idx].y, 
        hits[min_idx].z - hits[max_idx].z
    };

    return true;
}