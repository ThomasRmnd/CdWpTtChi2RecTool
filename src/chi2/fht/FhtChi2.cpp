#include "chi2/fht/FhtChi2.hpp"

#include <numeric>

double FhtChi2::calculate(const_iterator first, const_iterator last, theo_const_iterator theo) {
    return std::transform_reduce(first, last, theo, 0.0, std::plus<double>(), [](const RecPmtProp& pmt, double value) {
        return std::pow((pmt.fht - value) * pmt.inv_res, 2.0);
    });
}