#include "strategy/Factory.hpp"

StrategyFactory::StrategyFactory(StrategyRegistry& registry) :
    m_registry(registry)
{}

DetectorType StrategyFactory::getDetectorType(const RecPmtTable& table) const {
    DetectorType dtype = DetectorType::NONE;
    for (RecPmtTable::const_iterator it = table.begin(); it != table.end(); ++it) {
        if (!it->used) continue;
        if ( (it->type & RecPmtType::PMT_CD) == it->type ) dtype |= DetectorType::CD;
        else if ( (it->type & RecPmtType::PMT_WP) == it->type ) dtype |= DetectorType::WP;
        else if ( (it->type & RecPmtType::PMT_TT) == it->type ) dtype |= DetectorType::TT;
    }
    return dtype;
}

ManualStrategyFactory::ManualStrategyFactory(StrategyRegistry& registry, ParamsType ptype) :
    StrategyFactory(registry),
    c_ptype(ptype)
{}

std::shared_ptr<Strategy> ManualStrategyFactory::construct(const RecPmtTable& table) {
    return m_registry.get({c_ptype, getDetectorType(table)});
}

AutomaticStrategyFactory::AutomaticStrategyFactory(StrategyRegistry& registry, const std::shared_ptr<Classifier>& classifier) :
    StrategyFactory(registry),
    m_classifier(classifier)
{
    // TODO: remove comment when Classifier is fully implemented
    // if (!m_classifier) throw std::runtime_error("AutomaticStrategyFactory::AutomaticStrategyFactory: [ERROR] Classifier is nullptr\n");
}

std::shared_ptr<Strategy> AutomaticStrategyFactory::construct(const RecPmtTable& table) {
    return m_registry.get({m_classifier->classify(table), getDetectorType(table)});
}