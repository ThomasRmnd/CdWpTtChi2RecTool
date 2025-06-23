#ifndef CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_

#include "chi2/Chi2.hpp"

#include "utils/Chi2.hpp"

// Compute raw $\chi^2$ for the TT method
class TtChi2 : public Chi2<TtMethodTag> {

public:

    TtChi2(double res_) : m_ires2{1 / (res_ * res_)} {}

    ~TtChi2() override = default;

    double calculate(const_iterator first, const_iterator last, theo_const_iterator theo) {
        return m_ires2 * chi2TtMethod(first, last, theo);
    }

private:

    double m_ires2;

};

#endif // CDWPTTCHI2RECTOOL_CHI2_TT_TTCHI2_HPP_