#ifndef CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_

#include "chi2/Chi2.hpp"

/**
 * @class TtChi2
 * 
 * @brief Derived class for raw \f$ \chi^2 \f$ calculation for the TT method.
 */
class TtChi2 : public Chi2<TtMethodTag> {

public:

    TtChi2(const std::string& name, double res);

    ~TtChi2() override = default;

    void configure(const SniperJSON& config) override;

    /**
     * @brief Calculate the raw \f$ \chi^2 \f$.
     * 
     * @param first The beginning iterator of the experimental data.
     * @param last The end iterator of the experimental data.
     * @param theo The beginning iterator of the expected data.
     * 
     * @return The raw \f$ \chi^2 \f$.
     * 
     * \f$ \chi^2 = \sum_{k} \left( \frac{ \lVert(\overrightarrow{p_{k, meas}} - \overrightarrow{p_{k, theo}}) \rVert }{\sigma} \right)^2 \f$
     * 
     * - Where:
     * 
     * - \f$ \overrightarrow{p_{k, meas}} \f$ is the experimental position of the \f$ k^{th} \f$ hit.
     * 
     * - \f$ \overrightarrow{p_{k, meas}} \f$ is the expected position of the \f$ k^{th} \f$ hits.
     * 
     * - \f$ \sigma \f$ is the spatial resolution.
     * 
     * The range beginning at `theo` must contain at least `last - first` elements.
     */
    double calculate(const_iterator first, const_iterator last, theo_const_iterator theo);

private:

    double m_ires2;

};

#endif // CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_