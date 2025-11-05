#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_CORRECTIONMAPLOOPESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_CORRECTIONMAPLOOPESTIMATOR_HPP_

#include <memory>

#include "estimator/Estimator.hpp"
#include "estimator/fht/map/CorrectionMap.hpp"

class CorrectionMapLoopEstimator : public Estimator {

public:

    CorrectionMapLoopEstimator(const std::string& name, const std::shared_ptr<Estimator>& esti, const std::vector<std::shared_ptr<CorrectionMap>>& maps, std::size_t nb_loop);

    virtual ~CorrectionMapLoopEstimator() override = default;

    void configure(const SniperJSON& config) override;
    bool initialize() override;
    bool finalize() override;

    ParamsType getIParamsType() override;
    ParamsType getOParamsType() override;

    bool estimate(RecPmtTable& table) override;

protected:

    std::shared_ptr<Estimator> m_esti;
    std::vector<std::shared_ptr<CorrectionMap>> m_maps;
    std::size_t m_nb_loop;

    std::vector<double> m_fhts;

    void copyTable(const RecPmtTable& table);
    bool applyCorrMap(RecPmtTable& table);
    void resetFht(RecPmtTable& table);

    friend class FhtTtCorrMapEstimator;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_CORRECTIONMAPLOOPESTIMATOR_HPP_