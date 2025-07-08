#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_TT_FHTTTCORRMAPESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_TT_FHTTTCORRMAPESTIMATOR_HPP_

#include "estimator/MinimizerEstimator.hpp"

#include <algorithm>
#include <set>

#include "estimator/fht/CorrectionMapEstimator.hpp"
#include "estimator/tt/TtMinimizerEstimator.hpp"

template<typename Params>
class FhtTtCorrMapEstimator : public MinimizerEstimator<FhtTtMethodTag> {

    struct ScoredHits {

        std::vector<vec3> hits;
        double score;
        std::vector<double> fvars;

        friend bool operator<(const ScoredHits& lhs, const ScoredHits& rhs) {
            return lhs.score < rhs.score;
        };

    };

public:

    FhtTtCorrMapEstimator(const std::string& name) :
        MinimizerEstimator<FhtTtMethodTag>(name),
        c_max_nb_hits(0), c_nb_loop(0)
    {};

    FhtTtCorrMapEstimator(const std::string& name, const std::shared_ptr<Optimizer>& opti, const std::shared_ptr<CostFunction<FhtTtMethodTag>>& func, const std::shared_ptr<CorrectionMapLoopEstimator<Params>>& fht_esti, const std::shared_ptr<TtMinimizerEstimator>& tt_esti, std::size_t max_nb_trks, std::size_t nb_loop) :
        MinimizerEstimator<FhtTtMethodTag>(name, opti, func),
        m_fht_esti(fht_esti),
        m_tt_esti(tt_esti),
        c_max_nb_hits(max_nb_trks),
        c_nb_loop(nb_loop),
        m_curr_nb_hits(0)
    {};

    virtual ~FhtTtCorrMapEstimator() override = default;

    bool initialize() override {
        if (c_max_nb_hits == 0) LogWarn << "The maximum number of TT hits is set to 0\n";
        if (c_nb_loop == 0) LogWarn << "The number of loops is set to 0\n";
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
    };

    bool finalize() override {
        if (!m_fht_esti) return false;
        if (!m_fht_esti->finalize()) return false;
        if (!m_tt_esti) return false;
        if (!m_tt_esti->finalize()) return false;
        return MinimizerEstimator<FhtTtMethodTag>::finalize();
    };

    ParamsType getIParamsType() override {
        if (this->m_func->getIParamsType() != m_fht_esti->getIParamsType()) {
            LogError << "Cost function and correction maps estimator have different parameters type\n";
            return ParamsType::Unknown;
        }
        if (this->m_func->getIParamsType() != m_tt_esti->getIParamsType()) {
            LogError << "Cost function and TT estimator have different parameters type\n";
            return ParamsType::Unknown;
        }
        return this->m_func->getIParamsType();
    };

    ParamsType getOParamsType() override {
        if (this->m_func->getIParamsType() != m_fht_esti->getOParamsType()) {
            LogError << "Cost function and correction maps estimator have different parameters type\n";
            return ParamsType::Unknown;
        }
        if (this->m_func->getIParamsType() != m_tt_esti->get0ParamsType()) {
            LogError << "Cost function and TT estimator have different parameters type\n";
            return ParamsType::Unknown;
        }
        return this->m_func->getIParamsType();
    };

    virtual bool estimate(RecPmtTable& table) override {
        if (m_params.empty()) {
            LogError << "No parameters were set as input. Abording\n";
            return false;
        }
        m_scored_hits.clear();

        if (!m_fht_esti->setParams(this->m_ivars, this->m_steps, this->m_names)) return false;
        if (!m_fht_esti->estimate(table)) return false;
        this->m_params = m_fht_esti->getParams();
        this->m_cost = m_fht_esti->getCost();


        if (!m_tt_esti->setParams({0.0, 0.0, 0.0, 0.0, 0.0, -1.0}, {1.0, 1.0, 0.0, 1.0, 1.0, 0.0}, {"x", "y", "z", "dx", "dy", "dz"})) return false;
        if (!m_tt_esti->estimate(table)) {
            LogWarn << "Failed to reconstruct with TT: TT will not be used for this event\n";
            return true;
        }
        m_curr_nb_hits = std::min(c_max_nb_hits, m_tt_esti->getNbHits());
        if (!m_curr_nb_hits) {
            LogWarn << "No hits in the TT: TT will not be used for this event\n";
            return true;
        }

        // First Minimization with 3 points or more ----------

        if (!m_fht_esti->applyCorrMap(this->m_params, table)) return false;
        this->m_func->set(table);
        this->m_opti->setParams(this->m_params, this->m_steps, this->m_names);

        for (std::size_t k = 0; k < m_curr_nb_hits; ++k) {
            this->m_func->set(m_tt_esti->getIthHits(k));
            printHits(m_tt_esti->getIthHits(k));
            this->m_opti->optimize(*(this->m_func));
            m_scored_hits.push_back({m_tt_esti->getIthHits(k), this->m_opti->getCost(), this->m_opti->getParams()});
        }

        m_fht_esti->resetFht(table);

        // Second Minimization with 3 points or more ----------

        for (std::size_t k = 1; k < c_nb_loop; ++k) {
            for (typename std::vector<ScoredHits>::iterator it = m_scored_hits.begin(); it != m_scored_hits.end(); ++it) {
                if (!m_fht_esti->applyCorrMap(it->fvars, table)) return false;
                this->m_func->set(table);
                this->m_opti->setParams(it->fvars, this->m_steps, this->m_names);
                this->m_func->set(it->hits);
                printHits(it->hits);
                this->m_opti->operator()(*(this->m_func));
                it->score = this->m_opti->getCost();
                it->fvars = this->m_opti->getParams();
                m_fht_esti->resetFht(table);
            }
        }


        // ----------

        std::sort(m_scored_hits.begin(), m_scored_hits.end());
        for (const ScoredHits& hit : m_scored_hits) {
            LogDebug << "Hits: ";
            for (const vec3& h : hit.hits) {
                std::cout << "(" << h.x << ", " << h.y << ", " << h.z << "), ";
            }
            std::cout << "Cost: " << hit.score << std::endl;
        }
        this->m_func->set(m_scored_hits.front().hits);
        RecPmtTable dummy_empty_table;
        this->m_func->set(dummy_empty_table); // since the table is empty, this will only calculate the tt cost
        double tt_cost = this->m_func->operator()(m_scored_hits.front().fvars.data());
        LogDebug << "CD cost: " << m_scored_hits.front().score - tt_cost << ", TT cost: " << tt_cost << std::endl;

        if ( ((m_scored_hits.front().score - tt_cost) - this->m_cost) / this->m_cost < 0.175 && tt_cost < 0.05) {
            this->m_cost = m_scored_hits.front().score;
            this->m_params = m_scored_hits.front().fvars;
            LogDebug << "Final parameters: ";
            for (std::size_t k = 0; k < this->m_dim; ++k) {
                std::cout << this->m_names[k] << " = " << this->m_params[k] << ", ";
            }
            std::cout << "cost = " << this->m_cost << std::endl;
            return true;
        }

        LogInfo << "The cost is too high: " << m_scored_hits.front().score << "; we will switch to 2 points reconstruction." << std::endl;

        // Final Minimization with 2 points ----------

        if (!m_fht_esti->applyCorrMap(this->m_params, table)) return false;
        this->m_func->set(table);
        this->m_opti->setParams(this->m_params, this->m_steps, this->m_names);

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
                    this->m_func->set(curr_hits);
                    printHits(curr_hits);
                    this->m_opti->operator()(*(this->m_func));
                    m_scored_hits.push_back({curr_hits, this->m_opti->getCost(), this->m_opti->getParams()});
                } while (std::next_permutation(hits_mask.begin(), hits_mask.end()));
            }
        }

        m_fht_esti->resetFht(table);

        // Second Minimization with 2 points ----------

        for (std::size_t k = 1; k < c_nb_loop; ++k) {
            for (std::size_t i = nb_3hits; i < m_scored_hits.size(); ++i) {
                if (!m_fht_esti->applyCorrMap(m_scored_hits[i].fvars, table)) return false;
                this->m_func->set(table);
                this->m_opti->setParams(m_scored_hits[i].fvars, this->m_steps, this->m_names);
                this->m_func->set(m_scored_hits[i].hits);
                printHits(m_scored_hits[i].hits);
                this->m_opti->operator()(*(this->m_func));
                m_scored_hits[i].score = this->m_opti->getCost();
                m_scored_hits[i].fvars = this->m_opti->getParams();
                m_fht_esti->resetFht(table);
            }
        }

        std::sort(m_scored_hits.begin(), m_scored_hits.end());
        for (const ScoredHits& hit : m_scored_hits) {
            LogDebug << "Hits: ";
            for (const vec3& h : hit.hits) {
                std::cout << "(" << h.x << ", " << h.y << ", " << h.z << "), ";
            }
            std::cout << "Cost: " << hit.score << std::endl;
        }
        this->m_func->set(m_scored_hits.front().hits);
        this->m_func->set(dummy_empty_table); // since the table is empty, this will only calculate the tt cost
        tt_cost = this->m_func->operator()(m_scored_hits.front().fvars.data());
        LogDebug << "CD cost: " << m_scored_hits.front().score - tt_cost << ", TT cost: " << tt_cost << std::endl;

        this->m_cost = m_scored_hits.front().score;
        this->m_params = m_scored_hits.front().fvars;

        LogDebug << "Final parameters: ";
        for (std::size_t k = 0; k < this->m_dim; ++k) {
            std::cout << this->m_names[k] << " = " << this->m_params[k] << ", ";
        }
        std::cout << "cost = " << this->m_cost << std::endl;
        return true;
    };

protected:

    std::shared_ptr<CorrectionMapLoopEstimator<Params>> m_fht_esti;
    std::shared_ptr<TtMinimizerEstimator> m_tt_esti;
    const std::size_t c_max_nb_hits;
    const std::size_t c_nb_loop;
    std::size_t m_curr_nb_hits;
    std::vector<ScoredHits> m_scored_hits;

    void printHits(const std::vector<vec3>& hits) {
        LogDebug << "Hits: ";
        for (const vec3& hit : hits) {
            std::cout << '(' << hit.x << ", " << hit.y << ", " << hit.z << "), ";
        }
        std::cout << '\n';
    };

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_TT_FHTTTCORRMAPESTIMATOR_HPP_