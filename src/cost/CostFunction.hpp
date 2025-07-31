#ifndef CDWPTTCHI2RECTOOL_COST_COSTFUNCTION_HPP_
#define CDWPTTCHI2RECTOOL_COST_COSTFUNCTION_HPP_

#include <memory>

#include "SniperKernel/SniperLog.h"

#include "chi2/Chi2.hpp"
#include "predictor/Predictor.hpp"
#include "utils/Configurable.hpp"
#include "utils/Method.hpp"
#include "utils/TrackParams.hpp"

/**
 * @class CostFunction
 * @tparam _MethodTag Method tag
 * 
 * @brief Base class for cost calculation
 */
template<typename _MethodTag>
class CostFunction : public IParamsHandler, public Configurable {

    static_assert(std::is_base_of<MethodTag, _MethodTag>::value, "Tag must be derived from MethodTag");

public:

    typedef typename MethodTraits<_MethodTag>::vector_type vector_type;
    typedef typename MethodTraits<_MethodTag>::const_iterator const_iterator;
    typedef typename MethodTraits<_MethodTag>::theo_vector_type theo_vector_type;

    CostFunction(const std::string& name, const std::shared_ptr<Predictor<_MethodTag>>& pred, const std::shared_ptr<Chi2<_MethodTag>>& sim) :
        Configurable{name},
        m_pred{pred},
        m_chi2{sim}
    {}

    virtual ~CostFunction() = default;

    ParamsType getIParamsType() override {
        return m_pred->getIParamsType();
    }

    void configure(const SniperJSON& config) override {
        if (!config.valid()) return;
        m_chi2->configure(config);
    }

    // Set the experimental data
    virtual void set(const vector_type& data) = 0;

    /**
     * @brief 1. Calculate the expected data based on the track parameters
     * 2. Calculate the cost by comparing with the experimental data provided 
     * 
     * @param params the track parameters
     * 
     * @return Cost
     */
    virtual double operator()(const double* params) = 0;

protected:

    std::shared_ptr<Predictor<_MethodTag>> m_pred;
    std::shared_ptr<Chi2<_MethodTag>> m_chi2;

    const_iterator m_first;
    const_iterator m_last;
    theo_vector_type m_theo;

};

/**
 * @class CostFunction<FhtTtMethodTag>
 * 
 * @brief Base class for cost calculation, specialization with FhtTtMethodTag
 */
template<>
class CostFunction<FhtTtMethodTag> : public IParamsHandler, public Configurable {

public:

    CostFunction(
        const std::string& name,
        const std::shared_ptr<Predictor<FhtMethodTag>>& pred_fht, const std::shared_ptr<Chi2<FhtMethodTag>>& chi2_fht,
        const std::shared_ptr<Predictor<TtMethodTag>>& pred_tt, const std::shared_ptr<Chi2<TtMethodTag>>& chi2_tt
    ) :
        Configurable{name},
        m_pred_fht{pred_fht},
        m_chi2_fht{chi2_fht},
        m_pred_tt{pred_tt},
        m_chi2_tt{chi2_tt}
    {}

    virtual ~CostFunction() = default;

    ParamsType getIParamsType() override {
        if (m_pred_fht->getIParamsType() != m_pred_tt->getIParamsType()) {
            LogError << "Fht and Tt methods have different parameter types\n";
            return ParamsType::Unknown;
        }
        return m_pred_fht->getIParamsType();
    }

    void configure(const SniperJSON& config) override {
        if (!config.valid()) return;
        m_chi2_fht->configure(config);
        m_chi2_tt->configure(config);
    }

    // Set the FHT experimental data
    virtual void set(const RecPmtTable& table) = 0;

    // Set the TT experimental data
    virtual void set(const std::vector<vec3>& hits) = 0;
    
    /**
     * @brief 1. Calculate the expected data based on the track parameters
     * 2. Calculate the cost by comparing with the experimental data set 
     * 
     * @param params track parameters
     * 
     * @return Cost
     */
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