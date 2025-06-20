#ifndef CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_

#include "utils/method.hpp"

// Base class for raw $\chi^2$ calculation
template<typename _Tag>
class Chi2 {

    static_assert(std::is_base_of<method_tag, _Tag>::value, "Tag must derive from method_tag");

public:

    typedef typename method_traits<_Tag>::const_iterator const_iterator;
    typedef typename method_traits<_Tag>::theo_const_iterator theo_const_iterator;

    virtual ~Chi2() = default;

    virtual double calculate(const_iterator first, const_iterator last, theo_const_iterator theo) = 0;

};

#endif // CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_