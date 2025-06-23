#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_INITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_INITIALIZER_HPP_

#include "SniperKernel/ToolBase.h"
#include "utils/TrackParams.hpp"

#include "utils/Method.hpp"

template<typename _Tag>
class Initializer : public ToolBase, public OParamsHandler {

    static_assert(std::is_base_of<MethodTag, _Tag>::value, "Tag must be derived from MethodTag");

public:

    typedef typename MethodTraits<_Tag>::vector_type vector_type;

    using ToolBase::ToolBase;

    virtual ~Initializer() = default;

    virtual bool operator()(const vector_type& data) = 0;

    const std::vector<double>& getParams() const {
        return m_params;
    };

protected:

    std::vector<double> m_params;

};

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_INITIALIZER_HPP_