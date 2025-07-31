#ifndef CDWPTTCHI2RECTOOL_COST_FHT_TT_FHTTTCOSTFUNCTION_HPP_
#define CDWPTTCHI2RECTOOL_COST_FHT_TT_FHTTTCOSTFUNCTION_HPP_

#include "cost/CostFunction.hpp"

/**
 * @class FhtTtCostFunction
 * 
 * @brief Derived class for cost calculation for joint FHT method and TT method
 */
class FhtTtCostFunction : public CostFunction<FhtTtMethodTag> {

public:

    using CostFunction<FhtTtMethodTag>::CostFunction;

    ~FhtTtCostFunction() override = default;

    /**
     * @brief Set the FHT experimental data
     * 
     * @param data FHT experimental data vector
     */
    void set(const RecPmtTable& table) override;

    /**
     * @brief Set the TT experimental data
     * 
     * @param data TT experimental data vector
     */
    void set(const std::vector<vec3>& hits) override;
    
    /**
     * @brief 1. Calculate the expected data based on the track parameters
     * 2. Calculate the cost by comparing with the experimental data set 
     * 
     * @param params track parameters
     * 
     * @return Cost
     */
    double operator()(const double* params) override;

};

/**
 * @class FhtTtReducedCostFunction
 * 
 * @brief Derived class for cost calculation for joint FHT method and TT method
 */
class FhtTtReducedCostFunction : public CostFunction<FhtTtMethodTag> {

public:

    using CostFunction<FhtTtMethodTag>::CostFunction;

    ~FhtTtReducedCostFunction() override = default;

    // Set the FHT experimental data
    void set(const RecPmtTable& table) override;

    // Set the TT experimental data
    void set(const std::vector<vec3>& hits) override;

    /**
     * @brief 1. Calculate the expected data based on the track parameters
     * 2. Calculate the cost by comparing with the experimental data set 
     * 
     * @param params track parameters
     * 
     * @return Cost
     */
    double operator()(const double* params) override;

};

#endif // CDWPTTCHI2RECTOOL_COST_FHT_TT_FHTTTCOSTFUNCTION_HPP_