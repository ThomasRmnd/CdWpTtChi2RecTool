#ifndef CDWPTTCHI2RECTOOL_COST_FHT_FHTCOSTFUNCTION_HPP_
#define CDWPTTCHI2RECTOOL_COST_FHT_FHTCOSTFUNCTION_HPP_

#include "cost/CostFunction.hpp"

/**
 * @class FhtCostFunction
 * 
 * @brief Derived class for cost calculation for FHT method
 */
class FhtCostFunction : public CostFunction<FhtMethodTag> {

public:

    using CostFunction<FhtMethodTag>::CostFunction;

    ~FhtCostFunction() override = default;

    // Set the experimental data
    void set(const vector_type& table) override;

    /**
     * @brief 1. Calculate the expected data based on the track parameters
     * 2. Calculate the cost by comparing with the experimental data provided 
     * 3. Return the cost divided by NDF
     * 
     * @param params track parameters
     * 
     * @return Cost / NDF
     */
    double operator()(const double* params) override;

};

#endif // CDWPTTCHI2RECTOOL_COST_FHT_FHTCOSTFUNCTION_HPP_