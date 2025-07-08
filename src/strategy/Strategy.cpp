#include "strategy/Strategy.hpp"

#include "cost/fht/FhtCostFunction.hpp"

#include "estimator/ParamsModifierEstimator.hpp"
#include "estimator/fht/CorrectionMapEstimator.hpp"
#include "estimator/fht/FhtMinimizerEstimator.hpp"

#include "initializer/fht/ClusterBundleInitializer.hpp"
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
    t_end = t_start + length / constants::c;
    if (length < 0.0) {
        LogWarn << "Negative length: " << length << '\n';
        std::swap(start, end);
        std::swap(t_start, t_end);
    }
    tracks->addTrk(
        TVector3(start.x, start.y, start.z),
        TVector3(end.x, end.y, end.z),
        t_start, t_end,
        totpe, cost, 0
    );
}

// ################################################################################################
// ==================================== CD Stopping Strategy ======================================
// ################################################################################################

void CdStoppingStrategy::create() {
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

void CdStoppingStrategy::prepare() {
    getDefaultParams<ParamsType::SingleAcrylic>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdStoppingStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<ParamsType::SingleAcrylic>::set(params, t_start, start, dir, length);
    end = start + dir * length;
    t_end = t_start + length / constants::c;
    if (length < 0.0) {
        LogWarn << "Negative length: " << length << '\n';
        std::swap(start, end);
        std::swap(t_start, t_end);
    }
    tracks->addTrk(
        TVector3(start.x, start.y, start.z),
        TVector3(end.x, end.y, end.z),
        t_start, t_end,
        totpe, cost, 0
    );
}

// ################################################################################################
// ===================================== CD Double Strategy =======================================
// ################################################################################################

void CdDoubleStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterBundleInitializer>(
        "ClusterBundleInitializer", 1000.0, 5000.0
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
        g_pred_fht_cd_no_refl_ls_double,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "FhtMinimizerEstimator", g_opti, cost
    );
    m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<ParamsConstrainerEstimator<ParamsType::DoubleAcrylic>>(
        "ParamsConstrainerEstimator",
        std::make_shared<CorrectionMapLoopEstimator<ParamsType::DoubleAcrylic>>(
            "CorrectionMapLoopEstimator", g_opti, cost, 
            std::vector<std::shared_ptr<CorrectionMap<ParamsType::DoubleAcrylic>>>{g_corr_map_nnvt_double, g_corr_map_hamamatsu_double, g_corr_map_3inch_double}, 2u
        ),
        std::array<bool, 8>{false, false, false, false, false, false, true, true}
    );
    m_pipe->addStep(esti_corrmap_loop);
}

void CdDoubleStrategy::prepare() {
    getDefaultParams<ParamsType::DoubleAcrylic>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdDoubleStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start_1, start_2, dir_1, dir_2, end_1, end_2;
    double t_start_1, t_start_2, length_1, length_2, t_end_1, t_end_2;
    TrackSetterHelper<ParamsType::DoubleAcrylic>::set(params, t_start_1, start_1, dir_1, length_1);
    TrackSetterHelper<ParamsType::DoubleAcrylic>::set(params, t_start_2, start_2, dir_2, length_2);
    end_1 = start_1 + dir_1 * length_1;
    end_2 = start_2 + dir_2 * length_2;
    t_end_1 = t_start_1 + length_1 / constants::c;
    t_end_2 + t_start_2 + length_2 / constants::c;
    if (length_1 < 0.0) {
        LogWarn << "Negative length: " << length_1 << '\n';
        std::swap(start_1, end_1);
        std::swap(t_start_1, t_end_1);
    }
    if (length_2 < 0.0) {
        LogWarn << "Negative length: " << length_2 << '\n';
        std::swap(start_2, end_2);
        std::swap(t_start_2, t_end_2);
    }
    tracks->addTrk(
        TVector3(start_1.x, start_1.y, start_1.z),
        TVector3(end_1.x, end_1.y, end_1.z),
        t_start_1, t_end_1,
        totpe, cost, 0
    );
    tracks->addTrk(
        TVector3(start_2.x, start_2.y, start_2.z),
        TVector3(end_2.x, end_2.y, end_2.z),
        t_start_2, t_end_2,
        totpe, cost, 0
    );
}

// ################################################################################################
// ========================================= TT Strategy ==========================================
// ################################################################################################

// ################################################################################################
// ======================================== CD WP Strategy ========================================
// ################################################################################################

// ################################################################################################
// ======================================== CD TT Strategy ========================================
// ################################################################################################

// ################################################################################################
// ===================================== = CD WP TT Strategy ======================================
// ################################################################################################

// ################################################################################################
// ==================================== CD Water Phase Strategy ===================================
// ################################################################################################

// ################################################################################################
// =================================== CD WP Water Phase Strategy =================================
// ################################################################################################

// ################################################################################################
// ================================== CD WP TT Water Phase Strategy ===============================
// ################################################################################################