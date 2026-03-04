#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_HYBRIDLSWATERFHTPREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_HYBRIDLSWATERFHTPREDICTOR_HPP_

#include <memory>

#include "predictor/Predictor.hpp"
#include "utils/TrackParams.hpp"

/**
 * @class HybridLSWaterFhtPredictor
 * 
 * @brief Derived class to calculate expected data for FHT method either for a SingleTrack event
 * if the track goes through the acrylic and SingleCd otherwise
 */
class HybridLSWaterFhtPredictor : public Predictor<FhtMethodTag>, public TrackSetter<SingleCdParamsTag> {

public:

    HybridLSWaterFhtPredictor(const std::shared_ptr<Predictor<FhtMethodTag>>& pred_acrylic, const std::shared_ptr<Predictor<FhtMethodTag>>& pred_water);

    ~HybridLSWaterFhtPredictor() override = default;

    ParamsType getIParamsType() override;

    /**
     * @brief Calculate the expected data
     * 
     * @param first the beginning iterator of the experimental data 
     * @param last the end iterator of the experimental data
     * @param theo the beginning iterator of the expected data
     * @param params track parameters
     */
    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override;

protected:

    std::shared_ptr<Predictor<FhtMethodTag>> m_pred_acrylic;
    std::shared_ptr<Predictor<FhtMethodTag>> m_pred_water;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_HYBRIDLSWATERFHTPREDICTOR_HPP_