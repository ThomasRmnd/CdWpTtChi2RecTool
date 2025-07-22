#ifndef CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_

#include "chi2/Chi2.hpp"

/**
 * @class TtChi2
 * @brief Derived class for raw \f$ \chi^2 \f$ calculation for the TT method
 */
class TtChi2 : public Chi2<TtMethodTag> {

public:

    TtChi2(const std::string& name);
    TtChi2(const std::string& name, double res);

    ~TtChi2() override = default;

    void configure(const SniperJSON& config) override;

    /**
     * @brief Calculate the raw \f$ \chi^2 \f$ for the TT method
     * 
     * @param first the beginning iterator of the experimental data 
     * @param last the end iterator of the experimental data
     * @param theo the beginning iterator of the expected data
     * 
     * @return Raw \f$ \chi^2 \f$
     */
    double calculate(const_iterator first, const_iterator last, theo_const_iterator theo);

private:

    double m_ires2;

};

#endif // CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_