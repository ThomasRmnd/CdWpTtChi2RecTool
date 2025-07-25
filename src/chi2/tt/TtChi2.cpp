#include "chi2/tt/TtChi2.hpp"

#include <numeric>

TtChi2::TtChi2(const std::string& name, double res) :
    Chi2<TtMethodTag>{name},
    m_ires2{1 / (res * res)}
{}

void TtChi2::configure(const SniperJSON& config) {
    double res = getConfigValue<double>("SpatialRes", config);
    m_ires2 = 1 / (res * res);
}

double TtChi2::calculate(const_iterator first, const_iterator last, theo_const_iterator theo) {
    return m_ires2 * std::transform_reduce(first, last, theo, 0.0, std::plus<double>(), [](const vec3& hit, const vec3& value) {
        return mag2(hit - value);
    });
}