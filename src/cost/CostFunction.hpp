#ifndef CDWPTTCHI2RECTOOL_COST_COSTFUNCTION_HPP_
#define CDWPTTCHI2RECTOOL_COST_COSTFUNCTION_HPP_

#include "utils/TrackParams.hpp"

#include <memory>

#include "SniperKernel/SniperLog.h"

#include "chi2/Chi2.hpp"
#include "predictor/Predictor.hpp"

template<typename _Tag>
class CostFunction : public IParamsHandler {

    static_assert(std::is_base_of<MethodTag, _Tag>::value, "Tag must be derived from MethodTag");

public:

    typedef typename MethodTraits<_Tag>::vector_type vector_type;
    typedef typename MethodTraits<_Tag>::const_iterator const_iterator;
    typedef typename MethodTraits<_Tag>::theo_vector_type theo_vector_type;

    CostFunction(const std::shared_ptr<Predictor<_Tag>>& pred, const std::shared_ptr<Chi2<_Tag>>& sim) :
        m_pred(pred),
        m_chi2(sim)
    {}

    virtual ~CostFunction() = default;

    ParamsType getIParamsType() override {
        return m_pred->getIParamsType();
    }

    virtual void set(const vector_type& data) = 0;

    virtual double operator()(const double* params) = 0;

protected:

    std::shared_ptr<Predictor<_Tag>> m_pred;
    std::shared_ptr<Chi2<_Tag>> m_chi2;

    const_iterator m_first;
    const_iterator m_last;
    theo_vector_type m_theo;

};

template<>
class CostFunction<FhtTtMethodTag> : public IParamsHandler {

public:

    CostFunction(
        const std::shared_ptr<Predictor<FhtMethodTag>>& pred_fht, const std::shared_ptr<Chi2<FhtMethodTag>>& chi2_fht,
        const std::shared_ptr<Predictor<TtMethodTag>>& pred_tt, const std::shared_ptr<Chi2<TtMethodTag>>& chi2_tt
    ) :
        m_pred_fht(pred_fht),
        m_chi2_fht(chi2_fht),
        m_pred_tt(pred_tt),
        m_chi2_tt(chi2_tt)
    {}

    virtual ~CostFunction() = default;

    ParamsType getIParamsType() override {
        if (m_pred_fht->getIParamsType() != m_pred_tt->getIParamsType()) {
            LogError << "Fht and Tt methods have different parameter types\n";
            return ParamsType::Unknown;
        }
        return m_pred_fht->getIParamsType();
    }

    virtual void set(const RecPmtTable& table) = 0;
    virtual void set(const std::vector<vec3>& hits) = 0;
    
    virtual double operator()(const double* params) = 0;

protected:

    std::shared_ptr<Predictor<FhtMethodTag>> m_pred_fht;
    std::shared_ptr<Chi2<FhtMethodTag>> m_chi2_fht;
    std::shared_ptr<Predictor<TtMethodTag>> m_pred_tt;
    std::shared_ptr<Chi2<TtMethodTag>> m_chi2_tt;

    RecPmtTable::const_iterator m_ftable;
    RecPmtTable::const_iterator m_ltable;
    std::vector<double> m_fht_theo;

    std::vector<vec3>::const_iterator m_fhits;
    std::vector<vec3>::const_iterator m_lhits;
    std::vector<vec3> m_hits_theo;

};

#endif // CDWPTTCHI2RECTOOL_COST_COSTFUNCTION_HPP_