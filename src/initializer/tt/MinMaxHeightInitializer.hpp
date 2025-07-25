#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_TT_MINMAXHEIGHTINITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_TT_MINMAXHEIGHTINITIALIZER_HPP_

#include "initializer/Initializer.hpp"

/**
 * @class MinMaxHeightInitializer
 * 
 * @brief Derived class to calculate track parameters initial guess for TT method
 */
class MinMaxHeightInitializer : public Initializer<TtMethodTag> {

public:

    using Initializer<TtMethodTag>::Initializer;
    
    ~MinMaxHeightInitializer() override = default;

    ParamsType getOParamsType() override;

    /**
     * @brief Calculate the track parameters initial guess
     * 
     * @param data experimental data vector
     * 
     * @return Boolean whether the calculation is successfull or not
     */
    bool initiate(const vector_type& data) override;

};

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_TT_MINMAXHEIGHTINITIALIZER_HPP_