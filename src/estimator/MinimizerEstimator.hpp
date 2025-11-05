#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_MINIMIZERESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_MINIMIZERESTIMATOR_HPP_

#include <memory>

#include "cost/CostFunction.hpp"
#include "estimator/Estimator.hpp"
#include "estimator/optimizer/Optimizer.hpp"
#include "utils/Method.hpp"

template<typename _MethodTag>
class MinimizerEstimator : public Estimator {

    static_assert(std::is_base_of<MethodTag, _MethodTag>::value, "Tag must be derived from MethodTag");

public:

    MinimizerEstimator(const std::string& name) :
        Estimator(name),
        m_opti(nullptr),
        m_func(nullptr)
    {}

    MinimizerEstimator(const std::string& name, const std::shared_ptr<Optimizer>& opti, const std::shared_ptr<CostFunction<_MethodTag>>& func) :
        Estimator(name),
        m_opti(opti),
        m_func(func)
    {}

    virtual ~MinimizerEstimator() = default;

    virtual void configure(const SniperJSON& config) {
        if (!config.valid()) return;
        m_func->configure(config);
    }

    virtual bool initialize() override {
        if (!m_opti) {
            LogError << "Optimizer is not set\n";
            return false;
        }
        if (!m_func) {
            LogError << "Cost function is not set\n";
            return false;
        }
        if (getIParamsType() == ParamsType::Unknown) {
            LogError << "Problem with parameters compatibility, please check the cost function parameters type\n";
            return false;
        }
        return true;
    };

    virtual ParamsType getIParamsType() override {
        return m_func->getIParamsType();
    }

    virtual ParamsType getOParamsType() override {
        return m_func->getIParamsType();
    }

    virtual bool estimate(RecPmtTable& table) = 0;

protected:

    std::shared_ptr<Optimizer> m_opti;
    std::shared_ptr<CostFunction<_MethodTag>> m_func;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_MINIMIZERESTIMATOR_HPP_