#ifndef CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_
#define CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_

#include "utils/Configurable.hpp"
#include "utils/Method.hpp"

/**
 * @class Chi2
 * @tparam _MethodTag Method tag
 * 
 * @brief Base class for raw \f$ \chi^2 \f$ calculation.
 */
template<typename _MethodTag>
class Chi2 : public Configurable {

    static_assert(std::is_base_of<MethodTag, _MethodTag>::value, "Tag must derive from MethodTag");

public:

    typedef typename MethodTraits<_MethodTag>::const_iterator const_iterator;
    typedef typename MethodTraits<_MethodTag>::theo_const_iterator theo_const_iterator;

    using Configurable::Configurable;

    virtual ~Chi2() = default;

    /**
     * @brief Calculate the raw \f$ \chi^2 \f$.
     * 
     * @param first The beginning iterator of the experimental data.
     * @param last The end iterator of the experimental data.
     * @param theo The beginning iterator of the expected data.
     * 
     * @return The raw \f$ \chi^2 \f$.
     * 
     * The range beginning at `theo` must contain at least `last - first` elements.
     */
    virtual double calculate(const_iterator first, const_iterator last, theo_const_iterator theo) = 0;

};

#endif // CDWPTTCHI2RECTOOL_CHI2_CHI2_HPP_