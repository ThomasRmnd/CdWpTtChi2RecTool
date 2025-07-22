#include "chi2/fht/FhtChi2.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

#include "utils/Chi2.hpp"

DECLARE_TOOL(FhtChi2);

FhtChi2::FhtChi2(const std::string& name) :
    Chi2<FhtMethodTag>{name}
{}

double FhtChi2::calculate(const_iterator first, const_iterator last, theo_const_iterator theo) {
    return chi2FhtMethod(first, last, theo);
}