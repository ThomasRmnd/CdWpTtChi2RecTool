#include "chi2/tt/TtChi2.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

#include "utils/Chi2.hpp"

DECLARE_TOOL(TtChi2);

TtChi2::TtChi2(const std::string& name) :
    Chi2<TtMethodTag>{name}
{}

TtChi2::TtChi2(const std::string& name, double res) :
    Chi2<TtMethodTag>{name},
    m_ires2{1 / (res * res)}
{}

void TtChi2::configure(const SniperJSON& config) {
    double res = getConfigValue<double>("SpatialRes", config);
    m_ires2 = 1 / (res * res);
}

double TtChi2::calculate(const_iterator first, const_iterator last, theo_const_iterator theo) {
    return m_ires2 * chi2TtMethod(first, last, theo);
}