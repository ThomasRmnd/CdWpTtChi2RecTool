#include "cost/tt/TtCostFunction.hpp"

void TtCostFunction::set(const vector_type& hits) {
    m_first = hits.begin();
    m_last = hits.end();
    m_theo.resize(hits.size());
}

double TtCostFunction::operator()(const double* params) {
    m_pred->predict(m_first, m_last, m_theo.begin(), params);
    return m_chi2->calculate(m_first, m_last, m_theo.begin()) / static_cast<double>(3ul * std::distance(m_first, m_last) - g_track_type_to_size.at(getIParamsType()));
}