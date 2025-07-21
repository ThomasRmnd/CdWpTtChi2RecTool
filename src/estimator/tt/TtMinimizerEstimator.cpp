#include "estimator/tt/TtMinimizerEstimator.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(TtMinimizerEstimator);

TtMinimizerEstimator::TtMinimizerEstimator(const std::string& name) :
    MinimizerEstimator<TtMethodTag>(name, nullptr, nullptr),
    m_conv(nullptr),
    c_max_nb_hits(0ul),
    m_comb(nullptr),
    m_init(nullptr)
{}

TtMinimizerEstimator::TtMinimizerEstimator(const std::string& name, const std::shared_ptr<Optimizer>& opti, const std::shared_ptr<CostFunction<TtMethodTag>>& func, const std::shared_ptr<Converter>& conv, std::size_t max_nb_hits, const std::shared_ptr<Combinator>& comb, const std::shared_ptr<Initializer<TtMethodTag>>& init) :
    MinimizerEstimator<TtMethodTag>(name, opti, func),
    m_conv(conv),
    c_max_nb_hits(max_nb_hits),
    m_comb(comb),
    m_init(init)
{}

bool TtMinimizerEstimator::initialize() {
    if (!c_max_nb_hits) LogWarn << "The maximum number of TT hits is set to 0\n";
    if (!m_conv) {
        LogError << "Converter is not set\n";
        return false;
    }
    if (!m_conv->initialize()) return false;
    if (!m_comb) {
        LogError << "Combinator is not set\n";
        return false;
    }
    if (!m_comb->initialize()) return false;
    return MinimizerEstimator<TtMethodTag>::initialize();
}

bool TtMinimizerEstimator::finalize() {
    if (!m_comb) return false;
    if (!m_comb->finalize()) return false;
    if (!m_conv) return false;
    if (!m_conv->finalize()) return false;
    return MinimizerEstimator<TtMethodTag>::finalize();
}

ParamsType TtMinimizerEstimator::getIParamsType() {
    if (m_func->getIParamsType() != m_init->getOParamsType()) {
        LogError << "Cost function and initializer have different parameters type\n";
        return ParamsType::Unknown;
    }
    return m_func->getIParamsType();
}

ParamsType TtMinimizerEstimator::getOParamsType() {
    return getIParamsType();
}

std::size_t TtMinimizerEstimator::getNbHits() const { 
    return m_scored_hits.size(); 
}

const std::vector<vec3>& TtMinimizerEstimator::getIthHits(std::size_t i) const { 
    return m_scored_hits[i].hits; 
}

std::shared_ptr<CostFunction<TtMethodTag>>& TtMinimizerEstimator::getCostFunc() { 
    return m_func; 
}

bool TtMinimizerEstimator::estimate(RecPmtTable& table) {
    m_scored_hits.clear();
    
    RecPmtTable::const_iterator ftable = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return (pmt.type & (RecPmtType::PMT_CD | RecPmtType::PMT_WP)) == pmt.type; }).base();
    if (!m_conv->convert(ftable, table.end())) return false;
    if (c_max_nb_hits < m_conv->getHits().size()) {
        LogError << "The number of hits is too large: max -> " << c_max_nb_hits << ", actual -> " << m_conv->getHits().size() << '\n';
        return false;
    }
    if (!m_comb->combine(m_conv->getHits())) return false;
    
    std::vector<double> params;
    double cost = 0.;
    m_cost = std::numeric_limits<double>::max();

    for (const std::vector<vec3>& hits : *m_comb) {
        if (!m_init->initiate(hits)) return false;
        params = m_init->getParams();
        m_func->set(hits);
        if (!m_opti->setParams(params, m_steps, m_names)) return false;
        if (!m_opti->optimize(*m_func)) return false;
        cost = m_opti->getCost();
        if (cost < m_cost) {
            m_cost = cost;
            m_params = m_opti->getParams();
            m_opti->printParams();
        }
        m_scored_hits.push_back({hits, cost});
    }
    std::sort(m_scored_hits.begin(), m_scored_hits.end());
    
    return true;
};