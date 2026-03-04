#include "predictor/fht/HybridLSWaterFhtPredictor.hpp"

#include "utils/constants.hpp"

HybridLSWaterFhtPredictor::HybridLSWaterFhtPredictor(const std::shared_ptr<Predictor<FhtMethodTag>>& pred_acrylic, const std::shared_ptr<Predictor<FhtMethodTag>>& pred_water) :
    m_pred_acrylic(pred_acrylic),
    m_pred_water(pred_water)
{}

ParamsType HybridLSWaterFhtPredictor::getIParamsType() {
    if (m_pred_acrylic->getIParamsType() != ParamsType::SingleAcrylic) {
        LogError << "Acrylic predictor has not parameter type SingleAcrylic\n";
        return ParamsType::Unknown;
    }
    if (m_pred_water->getIParamsType() != ParamsType::SingleCd) {
        LogError << "Water Buffer predictor has not parameter type SingleCd\n";
        return ParamsType::Unknown;
    }
    return ParamsType::SingleCd;
}

void HybridLSWaterFhtPredictor::predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) {
    setTrack(params);
    double b_half = dot(m_dir, m_orig);
    double c = mag2(m_orig) - constants::r_acrylic * constants::r_acrylic;
    double disc_fourth = b_half * b_half - c;
    if (disc_fourth < 0.0) {
        m_pred_water->predict(first, last, theo, params); 
    } 
    else {
        double d = -b_half - std::sqrt(disc_fourth);
        vec3 orig_ = m_orig + d * m_dir;
        const double params_[5] = {params[0] + d / constants::c, theta(orig_), phi(orig_), params[3], params[4]};
        m_pred_acrylic->predict(first, last, theo, params_);
    }
}