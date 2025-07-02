#include "estimator/ParamsModifierEstimator.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(ParamsModifierEstimator);

bool ParamsModifierEstimator::estimate(RecPmtTable& table) {
    m_esti->setParams(m_params, m_steps, m_names);
    if (!m_esti->estimate(table)) return false;
    m_cost = m_esti->getCost();
    m_params = m_esti->getParams();
    return true;
}