#include "initializer/tt/MinMaxHeightInitializer.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(MinMaxHeightInitializer);

ParamsType MinMaxHeightInitializer::getOParamsType() {
    return ParamsType::SingleTt;
}

bool MinMaxHeightInitializer::initiate(const std::vector<vec3>& data) {
    std::size_t min_idx = 0, max_idx = 0;

    for (std::size_t k = 1; k < data.size(); ++k) {
        if (data[k].z < data[min_idx].z) min_idx = k;
        if (data[max_idx].z < data[k].z) max_idx = k;
    }

    if (min_idx == max_idx) {
        LogError << "Min and max hits are the same\n";
        return false;
    }

    m_params = std::vector<double>{
        data[max_idx].x, 
        data[max_idx].y, 
        data[max_idx].z,
        data[min_idx].x - data[max_idx].x, 
        data[min_idx].y - data[max_idx].y, 
        data[min_idx].z - data[max_idx].z
    };

    return true;
}