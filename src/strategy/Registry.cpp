#include "strategy/Registry.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(StrategyRegistry);

bool StrategyRegistry::initialize() {
    for (auto& [type, strat] : m_strats) {
        if (!strat) {
            LogError << "Strategy is not set. Abording\n";
            return false;
        }
        strat->create();
        std::shared_ptr<Pipeline> pipe = strat->pipeline();
        if (!pipe) {
            LogError << "Pipeline is not set for the strategy (" << static_cast<int>(strat->type.params) << ", " << static_cast<int>(strat->type.detector) << ")\n";
            return false;
        }
        if (!pipe->initialize()) return false;
    }
}

bool StrategyRegistry::finalize() {
    for (auto& [type, strat] : m_strats) {
        if (!strat) return false;
        std::shared_ptr<Pipeline> pipe = strat->pipeline();
        if (!pipe) return false;
        if (!pipe->finalize()) return false;
    }
}

void StrategyRegistry::book(const std::shared_ptr<Strategy>& strat) {
    if (m_strats.find(strat->type) != m_strats.end()) {
        LogWarn << "Already booked a strategy of type (" << static_cast<int>(strat->type.params) << ", " << static_cast<int>(strat->type.detector) << "). Overwritting\n";
    }
    m_strats[strat->type] = strat;
}

std::shared_ptr<Strategy> StrategyRegistry::get(const StrategyType& type) const {
    std::unordered_map<StrategyType, std::shared_ptr<Strategy>>::const_iterator it = m_strats.find(type);
    if (it == m_strats.end()) {
        std::cout << "StrategyRegistry::get: [WARN] No strategy found of type (" << static_cast<int>(type.params) << ", " << static_cast<int>(type.detector) << "). Skipping\n";
        return nullptr;
    }
    return it->second;
}