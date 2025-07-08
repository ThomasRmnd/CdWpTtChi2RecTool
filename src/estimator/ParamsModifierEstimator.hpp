#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_PARAMSMODIFIERESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_PARAMSMODIFIERESTIMATOR_HPP_

#include "estimator/Estimator.hpp"

#include <memory>

class ParamsModifierEstimator : public Estimator {

public:

    ParamsModifierEstimator(const std::string& name) : Estimator(name) {}
    ParamsModifierEstimator(const std::string& name, const std::shared_ptr<Estimator>& esti);

    virtual ~ParamsModifierEstimator() = default;

    virtual bool estimate(RecPmtTable& table);

protected:

    std::shared_ptr<Estimator> m_esti;

};

template<ParamsType _IPt>
class ParamsConstrainerEstimator : public ParamsModifierEstimator {

public:

    ParamsConstrainerEstimator(const std::string& name) : ParamsModifierEstimator(name) {}

    ParamsConstrainerEstimator(const std::string& name, const std::shared_ptr<Estimator> esti, const std::array<bool, ParamsTraits<_IPt>::size>& constrains) : 
        ParamsModifierEstimator(name, esti),
        m_constrains(constrains)
    {}

    ~ParamsConstrainerEstimator() override = default;

    bool estimate(RecPmtTable& table) override {
        if (m_params.size() != ParamsTraits<_IPt>::size) {
            LogError << "Parameters have wrong sizes, expected " << ParamsTraits<_IPt>::size << " but got " << m_params.size() << '\n';
            return false;
        }
        for (std::size_t i = 0; i < ParamsTraits<_IPt>::size; ++i) {
            if (m_constrains[i]) m_steps[i] = 0.0;
        }
        return ParamsModifierEstimator::estimate(table);
    }

    ParamsType get_IPtType() override {
        return _IPt;
    }

    ParamsType get_OPtType() override {
        return _IPt;
    }

private:

    std::array<bool, ParamsTraits<_IPt>::size> m_constrains;

};

namespace details {

void convert_params(std::vector<double>& params, std::vector<double>& steps, std::vector<std::string>& names, TrackSetterHelper<ParamsType::SingleAcrylic>, TrackSetterHelper<ParamsType::SingleStoppingAcrylic>) {
    params.push_back(1.0);
    steps.push_back(0.1);
    names.push_back("length");
}

} // namespace details

template<ParamsType _IPt, ParamsType _OPt>
class ParamsConverterEstimator : public ParamsModifierEstimator {

public:

    using ParamsModifierEstimator::ParamsModifierEstimator;

    ~ParamsConverterEstimator() override = default;

    bool estimate(RecPmtTable& table) override {
        if (m_params.size() != ParamsTraits<_IPt>::size) {
            LogError << "Parameters have wrong sizes, expected " << ParamsTraits<_IPt>::size << " but got " << m_params.size() << '\n';
            return false;
        }
        details::convert_params(m_params, m_steps, m_names, TrackSetterHelper<_IPt>{}, TrackSetterHelper<_OPt>{});
        return ParamsModifierEstimator::estimate(table);
    };

    ParamsType getIParamsType() override {
        return _IPt;
    };

    ParamsType getOParamsType() override {
        return _OPt;
    };

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_PARAMSMODIFIERESTIMATOR_HPP_