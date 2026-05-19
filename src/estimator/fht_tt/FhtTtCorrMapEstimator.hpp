#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_TT_FHTTTCORRMAPESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_TT_FHTTTCORRMAPESTIMATOR_HPP_

#include <algorithm>

#include "estimator/MinimizerEstimator.hpp"
#include "estimator/fht/CorrectionMapLoopEstimator.hpp"
#include "estimator/tt/TtMinimizerEstimator.hpp"

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

    FhtTtCorrMapEstimator(const std::string& name);
    FhtTtCorrMapEstimator(const std::string& name, const std::shared_ptr<Optimizer>& opti, const std::shared_ptr<CostFunction<FhtTtMethodTag>>& func, const std::shared_ptr<CorrectionMapLoopEstimator>& fht_esti, const std::shared_ptr<TtMinimizerEstimator>& tt_esti, std::size_t max_nb_trks, std::size_t nb_loop);

    ~FhtTtCorrMapEstimator() override = default;

    void configure(const SniperJSON& config) override;
    bool initialize() override;
    bool finalize() override;

    ParamsType getIParamsType() override;
    ParamsType getOParamsType() override;

    bool estimate(RecPmtTable& table) override;

protected:

    std::shared_ptr<CorrectionMapLoopEstimator> m_fht_esti;
    std::shared_ptr<TtMinimizerEstimator> m_tt_esti;
    std::size_t m_max_nb_hits;
    std::size_t m_nb_loop;
    
    std::size_t m_curr_nb_hits;
    std::vector<ScoredHits> m_scored_hits;

    void printHits(const std::vector<vec3>& hits) {
        std::ostringstream oss;
        oss << "Hits: ";
        for (const vec3& hit : hits) {
            oss << hit << ", ";
        }
        oss << '\n';
        LogDebug << oss.str();
    }

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_TT_FHTTTCORRMAPESTIMATOR_HPP_