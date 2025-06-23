#ifndef CDWPTTCHI2RECTOOL_COST_FHT_FHTCOSTFUNCTION_HPP_
#define CDWPTTCHI2RECTOOL_COST_FHT_FHTCOSTFUNCTION_HPP_

#include "cost/CostFunction.hpp"

class FhtCostFunction : public CostFunction<FhtMethodTag> {

public:

    using CostFunction<FhtMethodTag>::CostFunction;

    ~FhtCostFunction() override = default;

    void set(const vector_type& table) override {
        m_first = table.begin();
        m_last = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_CD | RecPmtType::PMT_WP); }).base();
        m_theo.resize(std::distance(m_first, m_last));
    }

    double operator()(const double* params) override {
        m_pred->predict(m_first, m_last, m_theo.begin(), params);
        return m_chi2->calculate(m_first, m_last, m_theo.begin()) / static_cast<double>(std::distance(m_first, m_last) - g_track_type_to_size.at(getIParamsType()));
    }

};

#endif // CDWPTTCHI2RECTOOL_COST_FHT_FHTCOSTFUNCTION_HPP_