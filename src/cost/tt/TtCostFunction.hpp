#ifndef CDWPTTCHI2RECTOOL_COST_TT_TTCOSTFUNCTION_HPP_
#define CDWPTTCHI2RECTOOL_COST_TT_TTCOSTFUNCTION_HPP_

#include "cost/CostFunction.hpp"

class TtCostFunction : public CostFunction<TtMethodTag> {

public:

    using CostFunction<TtMethodTag>::CostFunction;

    ~TtCostFunction() override = default;

    void set(const vector_type& hits) override {
        m_first = hits.begin();
        m_last = hits.end();
        m_theo.resize(hits.size());
    }

    double operator()(const double* params) override {
        m_pred->predict(m_first, m_last, m_theo.begin(), params);
        return m_chi2->calculate(m_first, m_last, m_theo.begin()) / static_cast<double>(3ul * std::distance(m_first, m_last) - g_track_type_to_size.at(getIParamsType()));
    }

};

#endif // CDWPTTCHI2RECTOOL_COST_TT_TTCOSTFUNCTION_HPP_