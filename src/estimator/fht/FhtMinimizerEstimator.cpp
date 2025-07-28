#include "estimator/fht/FhtMinimizerEstimator.hpp"

#include "SniperKernel/SniperLog.h"

bool FhtMinimizerEstimator::estimate(RecPmtTable& table) {
    if (m_params.empty()) {
        LogError << "No parameters were set as input. Abording\n";
        return false;
    }
    
    m_func->set(table);
    if (!m_opti->setParams(m_params, m_steps, m_names)) return false;
    if (!m_opti->optimize(*m_func)) return false;

    m_cost = m_opti->getCost();
    m_params = m_opti->getParams();

    m_opti->printParams();
    return true;
}