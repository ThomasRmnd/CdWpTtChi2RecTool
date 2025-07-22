#ifndef CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_

#include "utils/Configurable.hpp"
#include "utils/Method.hpp"

/**
 * @brief Base class for raw $\chi^2$ calculation
 */
template<typename _Tag>
class Chi2 : public Configurable {

    static_assert(std::is_base_of<MethodTag, _Tag>::value, "Tag must derive from MethodTag");

public:

    typedef typename MethodTraits<_Tag>::const_iterator const_iterator;
    typedef typename MethodTraits<_Tag>::theo_const_iterator theo_const_iterator;

    using Configurable::Configurable;

    virtual ~Chi2() = default;

    virtual double calculate(const_iterator first, const_iterator last, theo_const_iterator theo) = 0;

};

#endif // CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_