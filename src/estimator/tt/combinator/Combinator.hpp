#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_TT_COMBINATOR_COMBINATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_TT_COMBINATOR_COMBINATOR_HPP_

#include "SniperKernel/ToolBase.h"

#include "utils/vec3.hpp"

class Combinator : public ToolBase {

    typedef std::vector<std::vector<vec3>> PointCombination;

public:

    using ToolBase::ToolBase;

    virtual bool operator()(const std::vector<vec3>& hits) = 0;

    PointCombination::iterator begin() { return m_hits_comb.begin(); };
    PointCombination::const_iterator begin() const { return m_hits_comb.begin(); };
    PointCombination::const_iterator cbegin() const { return m_hits_comb.cbegin(); };
    PointCombination::iterator end() { return m_hits_comb.end(); };
    PointCombination::const_iterator end() const { return m_hits_comb.end(); };
    PointCombination::const_iterator cend() const { return m_hits_comb.cend(); };

protected:

    PointCombination m_hits_comb;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_TT_COMBINATOR_COMBINATOR_HPP_