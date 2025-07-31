#ifndef CDWPTTCHI2RECTOOL_COST_TT_TTCOSTFUNCTION_HPP_
#define CDWPTTCHI2RECTOOL_COST_TT_TTCOSTFUNCTION_HPP_

#include "cost/CostFunction.hpp"

/**
 * @class TtCostFunction
 * 
 * @brief Derived class for cost calculation for TT method
 */
class TtCostFunction : public CostFunction<TtMethodTag> {

public:

    using CostFunction<TtMethodTag>::CostFunction;

    ~TtCostFunction() override = default;

    // Set the experimental data
    void set(const vector_type& hits) override;

    /**
     * @brief 1. Calculate the expected data based on the track parameters
     * 2. Calculate the cost by comparing with the experimental data provided 
     * 
     * @param params track parameters
     * 
     * @return Cost
     */
    double operator()(const double* params) override;

};

#endif // CDWPTTCHI2RECTOOL_COST_TT_TTCOSTFUNCTION_HPP_