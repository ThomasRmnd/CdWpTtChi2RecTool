#include "estimator/fht_tt/FhtTtCorrMapEstimator.hpp"

FhtTtCorrMapEstimator::FhtTtCorrMapEstimator(const std::string& name) :
    MinimizerEstimator<FhtTtMethodTag>(name),
    m_max_nb_hits(0ul), m_nb_loop(0ul)
{}

FhtTtCorrMapEstimator::FhtTtCorrMapEstimator(const std::string& name, const std::shared_ptr<Optimizer>& opti, const std::shared_ptr<CostFunction<FhtTtMethodTag>>& func, const std::shared_ptr<CorrectionMapLoopEstimator>& fht_esti, const std::shared_ptr<TtMinimizerEstimator>& tt_esti, std::size_t max_nb_trks, std::size_t nb_loop) :
    MinimizerEstimator<FhtTtMethodTag>(name, opti, func),
    m_fht_esti(fht_esti),
    m_tt_esti(tt_esti),
    m_max_nb_hits(max_nb_trks),
    m_nb_loop(nb_loop),
    m_curr_nb_hits(0)
{}

void FhtTtCorrMapEstimator::configure(const SniperJSON& config) {
    if (!config.valid()) return;
    m_fht_esti->configure(config);
    m_tt_esti->configure(config);
    if (!setConfigValue(m_max_nb_hits, "MaxNumberHits", config)) return;
    if (!setConfigValue(m_nb_loop, "NumberLoops", config)) return;
}

bool FhtTtCorrMapEstimator::initialize() {
    if (m_max_nb_hits == 0) LogWarn << "The maximum number of TT hits is set to 0\n";
    if (m_nb_loop == 0) LogWarn << "The number of loops is set to 0\n";
    if (!m_fht_esti) {
        LogError << "FHT estimator is not set\n";
        return false;
    }
    if (!m_fht_esti->initialize()) return false;
    if (!m_tt_esti) {
        LogError << "TT estimator is not set\n";
        return false;
    }
    if (!m_tt_esti->initialize()) return false;
    return MinimizerEstimator<FhtTtMethodTag>::initialize();
}

bool FhtTtCorrMapEstimator::finalize() {
    if (!m_fht_esti) return false;
    if (!m_fht_esti->finalize()) return false;
    if (!m_tt_esti) return false;
    if (!m_tt_esti->finalize()) return false;
    return MinimizerEstimator<FhtTtMethodTag>::finalize();
}

ParamsType FhtTtCorrMapEstimator::getIParamsType() {
    if (m_func->getIParamsType() != m_fht_esti->getIParamsType()) {
        LogError << "Cost function and correction maps estimator have different parameters type\n";
        return ParamsType::Unknown;
    }
    return m_func->getIParamsType();
};

ParamsType FhtTtCorrMapEstimator::getOParamsType() {
    if (m_func->getIParamsType() != m_fht_esti->getOParamsType()) {
        LogError << "Cost function and correction maps estimator have different parameters type\n";
        return ParamsType::Unknown;
    }
    return m_func->getIParamsType();
};

bool FhtTtCorrMapEstimator::estimate(RecPmtTable& table) {
    if (m_params.empty()) {
        LogError << "No parameters were set as input. Abording\n";
        return false;
    }
    m_scored_hits.clear();

    if (!m_fht_esti->setParams(m_params, m_steps, m_names)) return false;
    if (!m_fht_esti->estimate(table)) return false;
    m_params = m_fht_esti->getParams();
    m_cost = m_fht_esti->getCost();


    if (!m_tt_esti->setParams({0.0, 0.0, 0.0, 0.0, 0.0, -1.0}, {1.0, 1.0, 0.0, 1.0, 1.0, 0.0}, {"x", "y", "z", "dx", "dy", "dz"})) return false;
    if (!m_tt_esti->estimate(table)) {
        LogWarn << "Failed to reconstruct with TT: TT will not be used for this event\n";
        return true;
    }
    m_curr_nb_hits = std::min(m_max_nb_hits, m_tt_esti->getNbHits());
    if (!m_curr_nb_hits) {
        LogWarn << "No hits in the TT: TT will not be used for this event\n";
        return true;
    }

    // First Minimization with 3 points or more ----------

    if (!m_fht_esti->applyCorrMap(table)) return false;
    m_func->set(table);
    m_opti->setParams(m_params, m_steps, m_names);

    for (std::size_t k = 0; k < m_curr_nb_hits; ++k) {
        m_func->set(m_tt_esti->getIthHits(k));
        printHits(m_tt_esti->getIthHits(k));
        m_opti->optimize(*m_func);
        m_scored_hits.push_back({m_tt_esti->getIthHits(k), m_opti->getCost(), m_opti->getParams()});
    }

    m_fht_esti->resetFht(table);

    // Second Minimization with 3 points or more ----------

    for (std::size_t k = 1; k < m_nb_loop; ++k) {
        for (typename std::vector<ScoredHits>::iterator it = m_scored_hits.begin(); it != m_scored_hits.end(); ++it) {
            if (!m_fht_esti->setParams(it->fvars, m_steps, m_names)) return false;
            if (!m_fht_esti->applyCorrMap(table)) return false;
            m_func->set(table);
            if (!m_opti->setParams(it->fvars, m_steps, m_names)) return false;
            m_func->set(it->hits);
            printHits(it->hits);
            m_opti->optimize(*m_func);
            it->score = m_opti->getCost();
            it->fvars = m_opti->getParams();
            m_fht_esti->resetFht(table);
        }
    }


    // ----------

    std::sort(m_scored_hits.begin(), m_scored_hits.end());
    for (const ScoredHits& hit : m_scored_hits) {
        printHits(hit.hits);
        std::cout << "Cost: " << hit.score << std::endl;
    }
    m_func->set(m_scored_hits.front().hits);
    RecPmtTable dummy_empty_table;
    m_func->set(dummy_empty_table); // since the table is empty, this will only calculate the tt cost
    double tt_cost = m_func->operator()(m_scored_hits.front().fvars.data());
    LogDebug << "CD cost: " << m_scored_hits.front().score - tt_cost << ", TT cost: " << tt_cost << std::endl;

    if ( ((m_scored_hits.front().score - tt_cost) - m_cost) / m_cost < 0.175 && tt_cost < 0.05) {
        m_cost = m_scored_hits.front().score;
        m_params = m_scored_hits.front().fvars;
        LogDebug << "Final parameters: ";
        for (std::size_t k = 0; k < m_size; ++k) {
            std::cout << m_names[k] << " = " << m_params[k] << ", ";
        }
        std::cout << "cost = " << m_cost << std::endl;
        return true;
    }

    LogInfo << "The cost is too high: " << m_scored_hits.front().score << "; we will switch to 2 points reconstruction." << std::endl;

    // Final Minimization with 2 points ----------

    if (!m_fht_esti->setParams(m_params, m_steps, m_names)) return false;
    if (!m_fht_esti->applyCorrMap(table)) return false;
    m_func->set(table);
    m_opti->setParams(m_params, m_steps, m_names);

    std::vector<bool> hits_mask;
    std::vector<vec3> curr_hits;

    std::size_t nb_3hits = m_scored_hits.size();
    for (std::size_t k = 0; k < nb_3hits; ++k) {
        printHits(m_scored_hits[k].hits);
        for (std::size_t i = m_scored_hits[k].hits.size() - 1; i > 1; --i) {
            hits_mask.assign(m_scored_hits[k].hits.size(), false);
            std::fill(hits_mask.end() - i, hits_mask.end(), true);
            do {
                LogDebug << "Mask: ";
                for (bool b : hits_mask) {
                    std::cout << b << ", ";
                }
                std::cout << std::endl;
                curr_hits.clear();
                for (std::size_t j = 0; j < hits_mask.size(); ++j) {
                    if (hits_mask[j]) curr_hits.push_back(m_scored_hits[k].hits[j]);
                }
                m_func->set(curr_hits);
                printHits(curr_hits);
                m_opti->optimize(*m_func);
                m_scored_hits.push_back({curr_hits, m_opti->getCost(), m_opti->getParams()});
            } while (std::next_permutation(hits_mask.begin(), hits_mask.end()));
        }
    }

    m_fht_esti->resetFht(table);

    // Second Minimization with 2 points ----------

    for (std::size_t k = 1; k < m_nb_loop; ++k) {
        for (std::size_t i = nb_3hits; i < m_scored_hits.size(); ++i) {
            if (!m_fht_esti->setParams(m_scored_hits[i].fvars, m_steps, m_names)) return false;
            if (!m_fht_esti->applyCorrMap(table)) return false;
            m_func->set(table);
            m_opti->setParams(m_scored_hits[i].fvars, m_steps, m_names);
            m_func->set(m_scored_hits[i].hits);
            printHits(m_scored_hits[i].hits);
            m_opti->optimize(*m_func);
            m_scored_hits[i].score = m_opti->getCost();
            m_scored_hits[i].fvars = m_opti->getParams();
            m_fht_esti->resetFht(table);
        }
    }

    std::sort(m_scored_hits.begin(), m_scored_hits.end());
    for (const ScoredHits& hit : m_scored_hits) {
        LogDebug << "Hits: ";
        for (const vec3& h : hit.hits) {
            std::cout << h << ", ";
        }
        std::cout << "Cost: " << hit.score << std::endl;
    }
    m_func->set(m_scored_hits.front().hits);
    m_func->set(dummy_empty_table); // since the table is empty, this will only calculate the tt cost
    tt_cost = m_func->operator()(m_scored_hits.front().fvars.data());
    LogDebug << "CD cost: " << m_scored_hits.front().score - tt_cost << ", TT cost: " << tt_cost << std::endl;

    m_cost = m_scored_hits.front().score;
    m_params = m_scored_hits.front().fvars;

    LogDebug << "Final parameters: ";
    for (std::size_t k = 0; k < m_size; ++k) {
        std::cout << m_names[k] << " = " << m_params[k] << ", ";
    }
    std::cout << "cost = " << m_cost << std::endl;
    return true;
}