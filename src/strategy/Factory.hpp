#ifndef CDWPTTCHI2RECTOOL_STRATEGY_FACTORY_HPP_
#define CDWPTTCHI2RECTOOL_STRATEGY_FACTORY_HPP_

#include "RecTools/IRecMuonTool.h"

#include "strategy/Registry.hpp"
#include "strategy/Strategy.hpp"
#include "strategy/classifier/Classifier.hpp"

class StrategyFactory {

public:

    StrategyFactory(const StrategyRegistry& registry);

    virtual ~StrategyFactory() = default;

    virtual std::shared_ptr<Strategy> construct(const RecPmtTable& table) = 0;

protected:

    const StrategyRegistry& c_registry;

};

class ManualStrategyFactory : public StrategyFactory {

public:

    ManualStrategyFactory(const StrategyRegistry& registry, ParamsType ptype);

    ~ManualStrategyFactory() override = default;

    std::shared_ptr<Strategy> construct(const RecPmtTable& table) override;

private:

    const ParamsType c_ptype;

};

class AutomaticStrategyFactory : public StrategyFactory {

public:

    AutomaticStrategyFactory(const StrategyRegistry& registry, const std::shared_ptr<Classifier>& classifier);

    ~AutomaticStrategyFactory() override = default;

    std::shared_ptr<Strategy> construct(const RecPmtTable& table) override;

private:

    std::shared_ptr<Classifier> m_classifier;

};

#endif // CDWPTTCHI2RECTOOL_STRATEGY_FACTORY_HPP_