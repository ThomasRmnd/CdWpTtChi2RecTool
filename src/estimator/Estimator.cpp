#include "estimator/Estimator.hpp"

#include "SniperKernel/SniperLog.h"

bool Estimator::setParams(const std::vector<double>& params, const std::vector<double>& steps, const std::vector<std::string>& names) {
    if (params.size() != steps.size() || params.size() != names.size()) {
        LogError << "Number of dimension from initial variables and steps is different\n";
        return false;
    }
    m_size = params.size();
    m_params = params;
    m_steps = steps;
    m_names = names;
    return true;
}

const std::vector<double>& Estimator::getParams() const {
    return m_params;
}

double Estimator::getCost() const {
    return m_cost;
}