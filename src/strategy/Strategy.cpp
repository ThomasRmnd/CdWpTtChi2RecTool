#include "strategy/Strategy.hpp"

#include "cost/fht/FhtCostFunction.hpp"

#include "estimator/fht/CorrectionMapEstimator.hpp"
#include "estimator/fht/FhtMinimizerEstimator.hpp"

#include "initializer/fht/ClusterMaxChargeInitializer.hpp"

#include "transformer/EarlyLateFhtTransformer.hpp"
#include "transformer/FhtChargeTholdTransformer.hpp"

#define DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(Base, name, Derived, ...) \
    std::shared_ptr<Base> g_##name##_single = std::make_shared<Derived<ParamsType::SingleAcrylic>>(__VA_ARGS__); \
    std::shared_ptr<Base> g_##name##_single_stopping = std::make_shared<Derived<ParamsType::SingleStoppingAcrylic>>(__VA_ARGS__); \
    std::shared_ptr<Base> g_##name##_double = std::make_shared<Derived<ParamsType::DoubleAcrylic>>(__VA_ARGS__);

#define DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(Base, name, Derived, ...) \
    std::shared_ptr<Base<ParamsType::SingleAcrylic>> g_##name##_single = std::make_shared<Derived<ParamsType::SingleAcrylic>>(__VA_ARGS__); \
    std::shared_ptr<Base<ParamsType::SingleStoppingAcrylic>> g_##name##_single_stopping = std::make_shared<Derived<ParamsType::SingleStoppingAcrylic>>(__VA_ARGS__); \
    std::shared_ptr<Base<ParamsType::DoubleAcrylic>> g_##name##_double = std::make_shared<Derived<ParamsType::DoubleAcrylic>>(__VA_ARGS__);

#include "chi2/fht/FhtChi2.hpp"
#include "chi2/tt/TtChi2.hpp"

std::shared_ptr<Chi2<FhtMethodTag>> g_chi2_fht = std::make_shared<FhtChi2>();
std::shared_ptr<Chi2<TtMethodTag>> g_chi2_tt = std::make_shared<TtChi2>(13.0);
std::shared_ptr<Chi2<TtMethodTag>> g_chi2_tt_joint = std::make_shared<TtChi2>(130.0);

#include "estimator/fht/map/CorrectionMap3d.hpp"
#include "estimator/fht/map/CorrParam.hpp"
#include "estimator/fht/map/WpTimeShiftCorrectionMap.hpp"

std::shared_ptr<CorrParam> g_corr_param_dist_proj_pmt_to_orig = std::make_shared<DistProjPmtToOrigCorrParam>();
std::shared_ptr<CorrParam> g_corr_param_angle = std::make_shared<AngleCorrParam>();
std::shared_ptr<CorrParam> g_corr_param_dist_track_to_center_squared = std::make_shared<DistTrackToCenterSquaredCorrParam>();

#include "estimator/optimizer/RootOptimizer.hpp"

std::shared_ptr<Optimizer> g_opti = std::make_shared<RootOptimizer>(1000000, 100000, 0.001);

// ################################################################################################
// ========================================= CD Strategy ==========================================
// ################################################################################################

void CdStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "ClusterMaxChargeInitializer", 5.0, 2.0, 5.0, 10.0, 9418.0, 0.9, 1.5, 10000.0
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
        "EarlyLateFhtTransformer", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    );
    m_pipe->addStep(trans_fht_spmt);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        g_pred_fht_cd_no_refl_ls_single,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "FhtMinimizerEstimator", g_opti, cost
    );
    m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<CorrectionMapLoopEstimator<ParamsType::SingleAcrylic>>(
        "CorrectionMapLoopEstimator", g_opti, cost, 
        std::vector<std::shared_ptr<CorrectionMap<ParamsType::SingleAcrylic>>>{g_corr_map_nnvt_single, g_corr_map_hamamatsu_single, g_corr_map_3inch_single}, 2u
    );
    m_pipe->addStep(esti_corrmap_loop);
}

void CdStrategy::prepare() {
    getDefaultParams<ParamsType::SingleAcrylic>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<ParamsType::SingleAcrylic>::set(params, t_start, start, dir, length);
    end = start + dir * length;
    t_end = t_start + length / 299.792458;
    tracks->addTrk(
        TVector3(start.x, start.y, start.z),
        TVector3(end.x, end.y, end.z),
        t_start, t_end,
        totpe, cost, 0 // TODO: Change this parameter to (cluster?)
    );
}