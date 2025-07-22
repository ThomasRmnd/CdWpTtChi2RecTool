#ifndef CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_

#include "SniperKernel/ToolBase.h"

#include "SniperKernel/SniperJSON.h"

class Configurable : public ToolBase {

public:

    using ToolBase::ToolBase;

    virtual ~Configurable() = default;

    virtual void configure(const SniperJSON& config) = 0;

};

#endif // CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_