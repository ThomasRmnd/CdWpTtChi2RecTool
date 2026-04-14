#include "cost/fht_tt/FhtTtCostFunction.hpp"

void FhtTtCostFunction::set(const RecPmtTable& table) {
    m_ftable = table.begin();
    m_ltable = std::find_if(
        table.rbegin(), table.rend(), 
        [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_CD | RecPmtType::PMT_WP); }
    ).base();
    m_fht_theo.resize(std::distance(m_ftable, m_ltable));
}

void FhtTtCostFunction::set(const std::vector<vec3>& hits) {
    m_fhits = hits.begin();
    m_lhits = hits.end();
    m_hits_theo.resize(hits.size());
}

double FhtTtCostFunction::operator()(const double* params) {
    m_pred_fht->predict(m_ftable, m_ltable, m_fht_theo.begin(), params);
    m_pred_tt->predict(m_fhits, m_lhits, m_hits_theo.begin(), params);
    return 
        ( m_chi2_fht->calculate(m_ftable, m_ltable, m_fht_theo.begin()) + m_chi2_tt->calculate(m_fhits, m_lhits, m_hits_theo.begin()) )
        / static_cast<double>(std::distance(m_ftable, m_ltable) + 3ul * std::distance(m_fhits, m_lhits) - g_track_type_to_size.at(getIParamsType()));
}

void FhtTtReducedCostFunction::set(const RecPmtTable& table) {
    m_ftable = table.begin();
    m_ltable = std::find_if(
        table.rbegin(), table.rend(), 
        [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_CD | RecPmtType::PMT_WP); }
    ).base();
    m_fht_theo.resize(std::distance(m_ftable, m_ltable));
}

void FhtTtReducedCostFunction::set(const std::vector<vec3>& hits) {
    m_fhits = hits.begin();
    m_lhits = hits.end();
    m_hits_theo.resize(hits.size());
}

double FhtTtReducedCostFunction::operator()(const double* params) {
    m_pred_fht->predict(m_ftable, m_ltable, m_fht_theo.begin(), params);
    m_pred_tt->predict(m_fhits, m_lhits, m_hits_theo.begin(), params);
    return 
        m_chi2_fht->calculate(m_ftable, m_ltable, m_fht_theo.begin()) 
        / static_cast<double>(std::distance(m_ftable, m_ltable) - g_track_type_to_size.at(getIParamsType())) 
        + m_chi2_tt->calculate(m_fhits, m_lhits, m_hits_theo.begin()) 
        / static_cast<double>(3ul * std::distance(m_fhits, m_lhits) - g_track_type_to_size.at(getIParamsType()));
}