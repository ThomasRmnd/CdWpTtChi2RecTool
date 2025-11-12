#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_CONVERTER_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_CONVERTER_HPP_

#include "utils/RecPmtProp.hpp"
#include "utils/Configurable.hpp"

class Converter : public Configurable {

public:

    using Configurable::Configurable;

    virtual ~Converter() = default;

    virtual bool convert(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) = 0;

    const std::vector<vec3>& hits() const {
        return m_hits;
    }

protected:

    std::vector<vec3> m_hits;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_CONVERTER_HPP_