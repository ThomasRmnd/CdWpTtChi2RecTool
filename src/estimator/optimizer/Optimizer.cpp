#include "estimator/optimizer/Optimizer.hpp"

#include "SniperKernel/SniperLog.h"

Optimizer::Optimizer(unsigned int max_calls, unsigned int max_iter, double tol) :
    c_max_calls(max_calls),
    c_max_iter(max_iter),
    c_tol(tol)
{}

bool Optimizer::setParams(const std::vector<double>& params, const std::vector<double>& steps, const std::vector<std::string>& names) {
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

const std::vector<double>& Optimizer::getParams() const {
    return m_params;
}

double Optimizer::getCost() const {
    return m_cost;
}

void Optimizer::printParams() {
    LogDebug << "Final parameters: ";
    for (std::size_t k = 0; k < m_size; ++k) {
        std::cout << m_names[k] << " = " << m_params[k] << ", ";
    }
    std::cout << "cost = " << m_cost << '\n';
}