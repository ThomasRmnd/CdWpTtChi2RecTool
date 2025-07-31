#ifndef CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_

#include "chi2/Chi2.hpp"

/**
 * @class FhtChi2
 * 
 * @brief Derived class for raw \f$ \chi^2 \f$ calculation for the FHT method.
 */
class FhtChi2 : public Chi2<FhtMethodTag> {

public:

    using Chi2<FhtMethodTag>::Chi2;

    ~FhtChi2() override = default;

    /**
     * @brief Calculate the raw \f$ \chi^2 \f$.
     * 
     * @param first The beginning iterator of the experimental data.
     * @param last The end iterator of the experimental data.
     * @param theo The beginning iterator of the expected data.
     * 
     * @return The raw \f$ \chi^2 \f$.
     * 
     * \f$ \chi^2 = \sum_{k} \left( \frac{t_{k, meas} - t_{k, theo}}{\sigma_{k}} \right)^{2} \f$
     * 
     * Where:
     * 
     * - \f$ t_{k, meas} \f$ is the experimental FHT of the \f$ k^{th} \f$ PMT.
     * 
     * - \f$ t_{k, theo} \f$ is the expected FHT of the \f$ k^{th} \f$ PMT.
     * 
     * - \f$ sigma_{k} \f$ is the time resolution of the \f$ k^{th} \f$ PMT.
     * 
     * 
     * The range beginning at `theo` must contain at least `last - first` elements.
     */
    double calculate(const_iterator first, const_iterator last, theo_const_iterator theo) override;

};

#endif // CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_
