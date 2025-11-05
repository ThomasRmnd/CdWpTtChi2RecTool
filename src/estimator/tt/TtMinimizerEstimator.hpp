#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_TT_TTMINIMIZERESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_TT_TTMINIMIZERESTIMATOR_HPP_

#include <cmath>

#include "estimator/MinimizerEstimator.hpp"
#include "estimator/tt/combinator/Combinator.hpp"
#include "estimator/tt/converter/Converter.hpp"
#include "initializer/Initializer.hpp"

class TtMinimizerEstimator : public MinimizerEstimator<TtMethodTag> {

    struct ScoredHits {

        std::vector<vec3> hits;
        double score;

        friend bool operator<(const ScoredHits& lhs, const ScoredHits& rhs) {
            return lhs.score < rhs.score;
        };

    };

public:

    TtMinimizerEstimator(const std::string& name, const std::shared_ptr<Optimizer>& opti, const std::shared_ptr<CostFunction<TtMethodTag>>& func, const std::shared_ptr<Converter>& conv, std::size_t max_nb_hits, const std::shared_ptr<Combinator>& comb, const std::shared_ptr<Initializer<TtMethodTag>>& init);

    ~TtMinimizerEstimator() override = default;

    void configure(const SniperJSON& config) override;
    bool initialize() override;
    bool finalize() override;

    ParamsType getIParamsType() override;
    ParamsType getOParamsType() override;

    std::size_t getNbHits() const;
    const std::vector<vec3>& getIthHits(std::size_t i) const;
    std::shared_ptr<CostFunction<TtMethodTag>>& getCostFunc();

    bool estimate(RecPmtTable& table) override;

protected:

    std::shared_ptr<Converter> m_conv;
    const std::size_t m_max_nb_hits;
    std::shared_ptr<Combinator> m_comb;
    std::shared_ptr<Initializer<TtMethodTag>> m_init;

    std::vector<ScoredHits> m_scored_hits;
    

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_TT_TTMINIMIZERESTIMATOR_HPP_