#include "cost/fht/FhtCostFunction.hpp"

void FhtCostFunction::set(const vector_type& table) {
    m_first = table.begin();
    m_last = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_CD | RecPmtType::PMT_WP); }).base();
    m_theo.resize(std::distance(m_first, m_last));
}

double FhtCostFunction::operator()(const double* params) {
    m_pred->predict(m_first, m_last, m_theo.begin(), params);
    return m_chi2->calculate(m_first, m_last, m_theo.begin()) / static_cast<double>(std::distance(m_first, m_last) - g_track_type_to_size.at(getIParamsType()));
}