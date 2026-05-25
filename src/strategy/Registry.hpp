#ifndef CDWPTTCHI2RECTOOL_STRATEGY_REGISTRY_HPP_
#define CDWPTTCHI2RECTOOL_STRATEGY_REGISTRY_HPP_

#include <unordered_map>

#include "strategy/Strategy.hpp"
#include "utils/Configurable.hpp"

class StrategyRegistry : public Configurable {

public:

    using Configurable::Configurable;

    void configure(const SniperJSON& config) override;
    bool initialize() override;
    bool finalize() override;

    void book(const std::shared_ptr<Strategy>& strat);
    std::shared_ptr<Strategy> get(const StrategyType& type);

private:

    SniperJSON m_config;
    std::unordered_map<StrategyType, std::shared_ptr<Strategy>> m_strats;

};

#endif // CDWPTTCHI2RECTOOL_STRATEGY_REGISTRY_HPP_