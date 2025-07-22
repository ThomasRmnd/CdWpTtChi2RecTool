#ifndef CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_

#include "chi2/Chi2.hpp"

/**
 * @class FhtChi2
 * @brief Derived class for raw \f$ \chi^2 \f$ calculation for the FHT method
 */
class FhtChi2 : public Chi2<FhtMethodTag> {

public:

    FhtChi2(const std::string& name);

    ~FhtChi2() override = default;

    /**
     * @brief Calculate the raw \f$ \chi^2 \f$ for the FHT method
     * 
     * @param first the beginning iterator of the experimental data 
     * @param last the end iterator of the experimental data
     * @param theo the beginning iterator of the expected data
     * 
     * @return Raw \f$ \chi^2 \f$
     */
    double calculate(const_iterator first, const_iterator last, theo_const_iterator theo) override;

};

#endif // CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_
