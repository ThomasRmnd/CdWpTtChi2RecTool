#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_CONVERTER_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_CONVERTER_HPP_

#include "SniperKernel/ToolBase.h"

#include "utils/RecPmtProp.hpp"

class Converter : public ToolBase {

public:

    using ToolBase::ToolBase;

    virtual ~Converter() = default;

    virtual bool convert(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) = 0;

    const std::vector<vec3>& getHits() {
        return m_hits;
    }

protected:

    std::vector<vec3> m_hits;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_CONVERTER_HPP_