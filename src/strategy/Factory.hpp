#ifndef CDWPTTCHI2RECTOOL_STRATEGY_FACTORY_HPP_
#define CDWPTTCHI2RECTOOL_STRATEGY_FACTORY_HPP_

#include "RecTools/IRecMuonTool.h"

#include "strategy/Registry.hpp"
#include "strategy/Strategy.hpp"
#include "strategy/classifier/Classifier.hpp"

class StrategyFactory {

public:

    StrategyFactory(StrategyRegistry& registry);

    virtual ~StrategyFactory() = default;

    virtual std::shared_ptr<Strategy> construct(const RecPmtTable& table) = 0;

protected:

    StrategyRegistry& m_registry;

    DetectorType getDetectorType(const RecPmtTable& table) const;

};

class ManualStrategyFactory : public StrategyFactory {

public:

    ManualStrategyFactory(StrategyRegistry& registry, ParamsType ptype);

    ~ManualStrategyFactory() override = default;

    std::shared_ptr<Strategy> construct(const RecPmtTable& table) override;

private:

    const ParamsType c_ptype;

};

class AutomaticStrategyFactory : public StrategyFactory {

public:

    AutomaticStrategyFactory(StrategyRegistry& registry, const std::shared_ptr<Classifier>& classifier);

    ~AutomaticStrategyFactory() override = default;

    std::shared_ptr<Strategy> construct(const RecPmtTable& table) override;

private:

    std::shared_ptr<Classifier> m_classifier;

};

#endif // CDWPTTCHI2RECTOOL_STRATEGY_FACTORY_HPP_