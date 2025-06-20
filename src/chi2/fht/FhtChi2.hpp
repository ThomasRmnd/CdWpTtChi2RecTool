#ifndef CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_

#include "chi2/Chi2.hpp"

#include "utils/chi2.hpp"

// Compute raw $\chi^2$ for the FHT method
class FhtChi2 : public Chi2<fht_method_tag> {

public:

    ~FhtChi2() override = default;

    double calculate(const_iterator first, const_iterator last, theo_const_iterator theo) override {
        return chi2_fht_method(first, last, theo);
    }

};

#endif // CDWPTTCHI2RECTOOL_CHI2_FHT_FHTCHI2_HPP_
