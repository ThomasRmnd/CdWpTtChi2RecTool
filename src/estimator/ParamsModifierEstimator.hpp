#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_PARAMSMODIFIERESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_PARAMSMODIFIERESTIMATOR_HPP_

#include "estimator/Estimator.hpp"

#include <memory>

class ParamsModifierEstimator : public Estimator {

public:

    ParamsModifierEstimator(const std::string& name, const std::shared_ptr<Estimator>& esti);

    virtual ~ParamsModifierEstimator() = default;

    virtual bool estimate(RecPmtTable& table);

protected:

    std::shared_ptr<Estimator> m_esti;

};

template<typename _ParamsTag>
class ParamsConstrainerEstimator : public ParamsModifierEstimator {

    static_assert(std::is_base_of<ParamsTag, _ParamsTag>::value, "Tag must derive from ParamsTag");

public:

    ParamsConstrainerEstimator(const std::string& name, const std::shared_ptr<Estimator> esti, const std::array<bool, ParamsTraits<_ParamsTag>::size>& constrains) : 
        ParamsModifierEstimator(name, esti),
        m_constrains(constrains)
    {}

    ~ParamsConstrainerEstimator() override = default;

    bool estimate(RecPmtTable& table) override {
        if (m_params.size() != ParamsTraits<_ParamsTag>::size) {
            LogError << "Parameters have wrong sizes, expected " << ParamsTraits<_ParamsTag>::size << " but got " << m_params.size() << '\n';
            return false;
        }
        for (std::size_t i = 0; i < ParamsTraits<_ParamsTag>::size; ++i) {
            if (m_constrains[i]) m_steps[i] = 0.0;
        }
        return ParamsModifierEstimator::estimate(table);
    }

    ParamsType getIParamsType() override {
        return ParamsTraits<_ParamsTag>::type;
    }

    ParamsType getOParamsType() override {
        return ParamsTraits<_ParamsTag>::type;
    }

private:

    std::array<bool, ParamsTraits<_ParamsTag>::size> m_constrains;

};

namespace details {

inline void convertParams(std::vector<double>& params, std::vector<double>& steps, std::vector<std::string>& names, SingleAcrylicParamsTag, SingleStoppingAcrylicParamsTag) {
    params.push_back(1.0);
    steps.push_back(0.1);
    names.push_back("length");
}

} // namespace details

template<typename _LhsParamsTag, typename _RhsParamsTag>
class ParamsConverterEstimator : public ParamsModifierEstimator {

    static_assert(std::is_base_of<ParamsTag, _LhsParamsTag>::value, "Tag must derive from ParamsTag");
    static_assert(std::is_base_of<ParamsTag, _RhsParamsTag>::value, "Tag must derive from ParamsTag");

public:

    using ParamsModifierEstimator::ParamsModifierEstimator;

    ~ParamsConverterEstimator() override = default;

    bool estimate(RecPmtTable& table) override {
        if (m_params.size() != ParamsTraits<_LhsParamsTag>::size) {
            LogError << "Parameters have wrong sizes, expected " << ParamsTraits<_LhsParamsTag>::size << " but got " << m_params.size() << '\n';
            return false;
        }
        details::convertParams(m_params, m_steps, m_names, _LhsParamsTag{}, _RhsParamsTag{});
        return ParamsModifierEstimator::estimate(table);
    };

    ParamsType getIParamsType() override {
        return ParamsTraits<_LhsParamsTag>::type;
    };

    ParamsType getOParamsType() override {
        return ParamsTraits<_RhsParamsTag>::type;
    };

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_PARAMSMODIFIERESTIMATOR_HPP_