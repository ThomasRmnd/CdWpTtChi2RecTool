#include "strategy/Strategy.hpp"

#include "chi2/fht/FhtChi2.hpp"
#include "chi2/tt/TtChi2.hpp"

#include "cost/fht/FhtCostFunction.hpp"
#include "cost/fht_tt/FhtTtCostFunction.hpp"
#include "cost/tt/TtCostFunction.hpp"

#include "estimator/ParamsModifierEstimator.hpp"
#include "estimator/fht/CorrectionMapEstimator.hpp"
#include "estimator/fht/FhtMinimizerEstimator.hpp"
#include "estimator/fht/map/CorrectionMap3d.hpp"
#include "estimator/fht/map/CorrParam.hpp"
#include "estimator/fht/map/WpTimeShiftCorrectionMap.hpp"
#include "estimator/fht_tt/FhtTtCorrMapEstimator.hpp"
#include "estimator/optimizer/RootOptimizer.hpp"
#include "estimator/tt/TtMinimizerEstimator.hpp"
#include "estimator/tt/combinator/MaskHeightCartesianProdCombinator.hpp"
#include "estimator/tt/converter/FuzeNeighborConverter.hpp"

#include "initializer/fht/ClusterBundleInitializer.hpp"
#include "initializer/fht/ClusterMaxChargeInitializer.hpp"
#include "initializer/fht/WaterPhaseInitializer.hpp"
#include "initializer/tt/MinMaxHeightInitializer.hpp"

#include "predictor/fht/CdFhtPredictor.hpp"
#include "predictor/fht/CdWpFhtPredictor.hpp"
#include "predictor/fht/WpFhtPredictor.hpp"
#include "predictor/fht/fht/CdFht.hpp"
#include "predictor/fht/fht/WpFht.hpp"
#include "predictor/tt/TtPredictor.hpp"

#include "transformer/CalibrationTransformer.hpp"
#include "transformer/EarlyLateFhtTransformer.hpp"
#include "transformer/FhtChargeTholdTransformer.hpp"
#include "transformer/TtCrossTalkTransformer.hpp"
#include "transformer/WaterPhaseTransformer.hpp"
#include "transformer/WpGeomTimeTransformer.hpp"

#define __USE_CCA_FILES_CORRECTION_MAP__

#ifdef __USE_CCA_FILES_CORRECTION_MAP__
    #define CORRECTION_MAP_FILENAME "~/J22.2.0-rc2/junosw/Examples/CorrectionMaps/share/CorrectionMap_CdWp_10010020_L2_5050_Shifted_2D_CCA.root"
    #define CORRECTION_MAP_MAPNAME_NNVT "NNVTPmtsCorrectionMap_FhtChargeThold20"
    #define CORRECTION_MAP_MAPNAME_HAMAMATSU "HamamatsuPmtsCorrectionMap_FhtChargeThold20"
    #define CORRECTION_MAP_MAPNAME_3INCH "3inchPmtsCorrectionMap_EarlyLateFht135"
    #define CORRECTION_MAP_MAPNAME_WP "WpPmtsCorrectionMap_TimeGeom"
#else
    #define CORRECTION_MAP_FILENAME "/junofs/users/traymond/data/CorrectionMaps/CorrectionMap_CdWp_10010020_L2_5050_Shifted_2D_CCA.root"
    #define CORRECTION_MAP_MAPNAME_NNVT "NNVTPmtsCorrectionMap_FhtChargeThold20"
    #define CORRECTION_MAP_MAPNAME_HAMAMATSU "HamamatsuPmtsCorrectionMap_FhtChargeThold20"
    #define CORRECTION_MAP_MAPNAME_3INCH "3inchPmtsCorrectionMap_EarlyLateFht135"
    #define CORRECTION_MAP_MAPNAME_WP "WpPmtsCorrectionMap_TimeGeom"
#endif // __USE_CCA_FILES_CORRECTION_MAP__

#define DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(Base, name, Derived, ...) \
    std::shared_ptr<Base> g_##name##_single = std::make_shared<Derived<ParamsType::SingleAcrylic>>(__VA_ARGS__); \
    std::shared_ptr<Base> g_##name##_single_stopping = std::make_shared<Derived<ParamsType::SingleStoppingAcrylic>>(__VA_ARGS__); \
    std::shared_ptr<Base> g_##name##_double = std::make_shared<Derived<ParamsType::DoubleAcrylic>>(__VA_ARGS__);

#define DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(Base, name, Derived, ...) \
    std::shared_ptr<Base<ParamsType::SingleAcrylic>> g_##name##_single = std::make_shared<Derived<ParamsType::SingleAcrylic>>(__VA_ARGS__); \
    std::shared_ptr<Base<ParamsType::SingleStoppingAcrylic>> g_##name##_single_stopping = std::make_shared<Derived<ParamsType::SingleStoppingAcrylic>>(__VA_ARGS__); \
    std::shared_ptr<Base<ParamsType::DoubleAcrylic>> g_##name##_double = std::make_shared<Derived<ParamsType::DoubleAcrylic>>(__VA_ARGS__);

#define DEFINIT_GLOBAL_PREDICTOR(Base, name, Derived, fhtname) \
    std::shared_ptr<Base> g_##name##_single = std::make_shared<Derived<ParamsType::SingleAcrylic>>(g_##fhtname##_single); \
    std::shared_ptr<Base> g_##name##_single_stopping = std::make_shared<Derived<ParamsType::SingleStoppingAcrylic>>(g_##fhtname##_single_stopping); \
    std::shared_ptr<Base> g_##name##_double = std::make_shared<Derived<ParamsType::DoubleAcrylic>>(g_##fhtname##_double);

#define DEFINIT_GLOBAL_PREDICTOR_CDWP(Base, name, Derived, cdfhtname, wpfhtname) \
    std::shared_ptr<Base> g_##name##_single = std::make_shared<Derived<ParamsType::SingleAcrylic>>(g_##cdfhtname##_single, g_##wpfhtname##_single); \
    std::shared_ptr<Base> g_##name##_single_stopping = std::make_shared<Derived<ParamsType::SingleStoppingAcrylic>>(g_##cdfhtname##_single_stopping, g_##wpfhtname##_single_stopping); \
    std::shared_ptr<Base> g_##name##_double = std::make_shared<Derived<ParamsType::DoubleAcrylic>>(g_##cdfhtname##_double, g_##wpfhtname##_double);

std::shared_ptr<Chi2<FhtMethodTag>> g_chi2_fht = std::make_shared<FhtChi2>();
std::shared_ptr<Chi2<TtMethodTag>> g_chi2_tt = std::make_shared<TtChi2>(13.0);
std::shared_ptr<Chi2<TtMethodTag>> g_chi2_tt_joint = std::make_shared<TtChi2>(130.0);

std::shared_ptr<CorrParam> g_corr_param_dist_proj_pmt_to_orig = std::make_shared<DistProjPmtToOrigCorrParam>();
std::shared_ptr<CorrParam> g_corr_param_angle = std::make_shared<AngleCorrParam>();
std::shared_ptr<CorrParam> g_corr_param_dist_track_to_center_squared = std::make_shared<DistTrackToCenterSquaredCorrParam>();

DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CorrectionMap, corr_map_nnvt, CorrectionMap3d, 
    "CorrectionMap3d", RecPmtType::PMT_20INCH_NNVT | RecPmtType::PMT_20INCH_HIGHQENNVT,
    CORRECTION_MAP_FILENAME, CORRECTION_MAP_MAPNAME_NNVT,
    g_corr_param_dist_proj_pmt_to_orig, g_corr_param_angle, g_corr_param_dist_track_to_center_squared
)
DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CorrectionMap, corr_map_hamamatsu, CorrectionMap3d,
    "CorrectionMap3d", RecPmtType::PMT_20INCH_HAMAMATSU,
    CORRECTION_MAP_FILENAME, CORRECTION_MAP_MAPNAME_HAMAMATSU,
    g_corr_param_dist_proj_pmt_to_orig, g_corr_param_angle, g_corr_param_dist_track_to_center_squared
)
DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CorrectionMap, corr_map_3inch, CorrectionMap3d,
    "CorrectionMap3d", RecPmtType::PMT_3INCH,
    CORRECTION_MAP_FILENAME, CORRECTION_MAP_MAPNAME_3INCH,
    g_corr_param_dist_proj_pmt_to_orig, g_corr_param_angle, g_corr_param_dist_track_to_center_squared
)
DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CorrectionMap, corr_map_wp, WpTimeShiftCorrectionMap,
    "WpTimeShiftCorrectionMap",
    CORRECTION_MAP_FILENAME, CORRECTION_MAP_MAPNAME_WP
)

std::shared_ptr<Optimizer> g_opti = std::make_shared<RootOptimizer>(1000000, 100000, 0.001);

DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CdFht, cd_fht, NoRefractionLsCdFht)
DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(WpFht, wp_fht, NoGeomWpFht)

DEFINIT_GLOBAL_PREDICTOR(Predictor<FhtMethodTag>, pred_fht_cd_no_refr_ls, CdFhtPredictor, cd_fht)
DEFINIT_GLOBAL_PREDICTOR(Predictor<FhtMethodTag>, pred_fht_wp_no_hit, WpFhtPredictor, wp_fht)
DEFINIT_GLOBAL_PREDICTOR_CDWP(Predictor<FhtMethodTag>, pred_fht_no_refl_ls_no_hit, CdWpFhtPredictor, cd_fht, wp_fht)
DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(Predictor<TtMethodTag>, pred_tt, TtPredictor)

std::shared_ptr<Transformer> g_trans_calib_hama = std::make_shared<CalibrationTransformer>("CalibrationTransformer", RecPmtType::PMT_20INCH_HAMAMATSU, 0.0);
std::shared_ptr<Transformer> g_trans_calib_nnvt = std::make_shared<CalibrationTransformer>("CalibrationTransformer", RecPmtType::PMT_20INCH_NNVT, 0.0);
std::shared_ptr<Transformer> g_trans_calib_highqe = std::make_shared<CalibrationTransformer>("CalibrationTransformer", RecPmtType::PMT_20INCH_HIGHQENNVT, 0.0);
std::shared_ptr<Transformer> g_trans_calib_spmt = std::make_shared<CalibrationTransformer>("CalibrationTransformer", RecPmtType::PMT_3INCH, 0.0);
std::shared_ptr<Transformer> g_trans_calib_wp = std::make_shared<CalibrationTransformer>("CalibrationTransformer", RecPmtType::PMT_WP, 0.0);

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
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

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
        g_pred_fht_cd_no_refr_ls_single,
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
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

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
        g_pred_fht_cd_no_refr_ls_single,
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
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

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
        g_pred_fht_cd_no_refr_ls_double,
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
    t_end_2 = t_start_2 + length_2 / constants::c;
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

void TtStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("Pipeline");

    // ======================================= Transformer ========================================
    std::shared_ptr<Transformer> trans_tt = std::make_shared<TtCrossTalkTransformer>(
        "TtCrossTalkTransformer"
    );
    m_pipe->addStep(trans_tt);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<Estimator> esti_tt = std::make_shared<TtMinimizerEstimator>(
        "TtMinimizerEstimator",
        g_opti,
        std::make_shared<TtCostFunction>(
            g_pred_tt_single_tt,
            g_chi2_tt
        ),
        std::make_shared<FuzeNeighborConverter>("FuzeNeighborConverter"), 25,
        std::make_shared<MaskHeightCartesianProdCombinator>("MaskHeightCartesianProdCombinator", 6),
        std::make_shared<MinMaxHeightInitializer>("MinMaxHeightInitializer")
    );
    m_pipe->addStep(esti_tt);
}

void TtStrategy::prepare() {
    getDefaultParams<ParamsType::SingleTt>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void TtStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<ParamsType::SingleTt>::set(params, t_start, start, dir, length);
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
// ======================================== CD WP Strategy ========================================
// ################################################################################################

void CdWpStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "ClusterMaxChargeInitializer", 5.0, 2.0, 5.0, 10.0, 9418.0, 0.9, 1.5, 10000.0
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);
    m_pipe->addStep(g_trans_calib_wp);

    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
        "EarlyLateFhtTransformer", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    );
    m_pipe->addStep(trans_fht_spmt);

    std::shared_ptr<Transformer> trans_q_wp = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_WP, 10.0
    );
    m_pipe->addStep(trans_q_wp);

    std::shared_ptr<Transformer> trans_geomtime_wp = std::make_shared<WpGeomTimeTransformer>(
        "WpGeomTimeTransformer", 50.0, 100.0, 0.4, 0.8, 10.0, 0.05
    );
    m_pipe->addStep(trans_geomtime_wp);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        g_pred_fht_no_refr_ls_no_hit_single,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "FhtMinimizerEstimator", g_opti, cost
    );
    m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<CorrectionMapLoopEstimator<ParamsType::SingleAcrylic>>(
        "CorrectionMapLoopEstimator", g_opti, cost, 
        std::vector<std::shared_ptr<CorrectionMap<ParamsType::SingleAcrylic>>>{g_corr_map_nnvt_single, g_corr_map_hamamatsu_single, g_corr_map_3inch_single, g_corr_map_wp_single}, 2u
    );
    m_pipe->addStep(esti_corrmap_loop);
}

void CdWpStrategy::prepare() {
    getDefaultParams<ParamsType::SingleAcrylic>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdWpStrategy::save(RecTrks* tracks, double totpe) {
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
// ======================================== CD TT Strategy ========================================
// ################################################################################################

void CdTtStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "ClusterMaxChargeInitializer", 5.0, 2.0, 5.0, 10.0, 9418.0, 0.9, 1.5, 10000.0
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
        "EarlyLateFhtTransformer", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    );
    m_pipe->addStep(trans_fht_spmt);

    std::shared_ptr<Transformer> trans_tt = std::make_shared<TtCrossTalkTransformer>(
        "TtCrossTalkTransformer"
    );
    m_pipe->addStep(trans_tt);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        g_pred_fht_cd_no_refr_ls_single,
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

    std::shared_ptr<CostFunction<TtMethodTag>> tt_cost = std::make_shared<TtCostFunction>(
        std::make_shared<TtPredictor<ParamsType::SingleTt>>(),
        g_chi2_tt
    );

    std::shared_ptr<Estimator> esti_tt = std::make_shared<TtMinimizerEstimator>(
        "TtMinimizerEstimator",
        g_opti, tt_cost,
        std::make_shared<FuzeNeighborConverter>("FuzeNeighborConverter"), 25,
        std::make_shared<MaskHeightCartesianProdCombinator>("MaskHeightCartesianProdCombinator", 6),
        std::make_shared<MinMaxHeightInitializer>("MinMaxHeightInitializer")
    );

    std::shared_ptr<CostFunction<FhtTtMethodTag>> fht_tt_cost = std::make_shared<FhtTtCostFunction>(
        g_pred_fht_cd_no_refr_ls_single,
        g_chi2_fht,
        g_pred_tt_single,
        g_chi2_tt_joint
    );

    std::shared_ptr<Estimator> esti_fhttt_corrmap = std::make_shared<FhtTtCorrMapEstimator<ParamsType::SingleAcrylic>>(
        "FhtTtCorrectionMapEstimator",
        g_opti, fht_tt_cost, 
        std::dynamic_pointer_cast<CorrectionMapLoopEstimator<ParamsType::SingleAcrylic>>(esti_corrmap_loop), 
        std::dynamic_pointer_cast<TtMinimizerEstimator>(esti_tt),
        10, 2
    );
    m_pipe->addStep(esti_fhttt_corrmap);
}

void CdTtStrategy::prepare() {
    getDefaultParams<ParamsType::SingleAcrylic>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdTtStrategy::save(RecTrks* tracks, double totpe) {
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
// ===================================== = CD WP TT Strategy ======================================
// ################################################################################################

void CdWpTtStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "ClusterMaxChargeInitializer", 5.0, 2.0, 5.0, 10.0, 9418.0, 0.9, 1.5, 10000.0
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);
    m_pipe->addStep(g_trans_calib_wp);

    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
        "EarlyLateFhtTransformer", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    );
    m_pipe->addStep(trans_fht_spmt);

    std::shared_ptr<Transformer> trans_q_wp = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_WP, 10.0
    );
    m_pipe->addStep(trans_q_wp);

    std::shared_ptr<Transformer> trans_geomtime_wp = std::make_shared<WpGeomTimeTransformer>(
        "WpGeomTimeTransformer", 50.0, 100.0, 0.4, 0.8, 10.0, 0.05
    );
    m_pipe->addStep(trans_geomtime_wp);

    std::shared_ptr<Transformer> trans_tt = std::make_shared<TtCrossTalkTransformer>(
        "TtCrossTalkTransformer"
    );
    m_pipe->addStep(trans_tt);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        g_pred_fht_no_refr_ls_no_hit_single,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "FhtMinimizerEstimator", g_opti, cost
    );
    m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<CorrectionMapLoopEstimator<ParamsType::SingleAcrylic>>(
        "CorrectionMapLoopEstimator", g_opti, cost, 
        std::vector<std::shared_ptr<CorrectionMap<ParamsType::SingleAcrylic>>>{g_corr_map_nnvt_single, g_corr_map_hamamatsu_single, g_corr_map_3inch_single, g_corr_map_wp_single}, 2
    );

    std::shared_ptr<CostFunction<TtMethodTag>> tt_cost = std::make_shared<TtCostFunction>(
        g_pred_tt_single_tt,
        g_chi2_tt
    );

    std::shared_ptr<Estimator> esti_tt = std::make_shared<TtMinimizerEstimator>(
        "TtMinimizerEstimator",
        g_opti, tt_cost,
        std::make_shared<FuzeNeighborConverter>("FuzeNeighborConverter"), 25,
        std::make_shared<MaskHeightCartesianProdCombinator>("MaskHeightCartesianProdCombinator", 6),
        std::make_shared<MinMaxHeightInitializer>("MinMaxHeightInitializer")
    );

    std::shared_ptr<CostFunction<FhtTtMethodTag>> fht_tt_cost = std::make_shared<FhtTtCostFunction>(
        g_pred_fht_no_refr_ls_no_hit_single,
        g_chi2_fht,
        g_pred_tt_single,
        g_chi2_tt_joint
    );

    std::shared_ptr<Estimator> esti_fhttt_corrmap = std::make_shared<FhtTtCorrMapEstimator<ParamsType::SingleAcrylic>>(
        "FhtTtCorrectionMapEstimator",
        g_opti, fht_tt_cost, 
        std::dynamic_pointer_cast<CorrectionMapLoopEstimator<ParamsType::SingleAcrylic>>(esti_corrmap_loop), 
        std::dynamic_pointer_cast<TtMinimizerEstimator>(esti_tt),
        10, 2
    );
    m_pipe->addStep(esti_fhttt_corrmap);
}

void CdWpTtStrategy::prepare() {
    getDefaultParams<ParamsType::SingleAcrylic>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdWpTtStrategy::save(RecTrks* tracks, double totpe) {
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
// ==================================== CD Water Phase Strategy ===================================
// ################################################################################################

void CdWaterPhaseStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<WaterPhaseInitializer>(
        "WaterPhaseInitializer", 200.0, 20.0, 0.25
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

    std::shared_ptr<Transformer> trans_hama_calib = std::make_shared<CalibrationTransformer>(
        "CalibrationTransformer", RecPmtType::PMT_20INCH_HAMAMATSU, -7.0
    );
    m_pipe->addStep(trans_hama_calib);
    
    std::shared_ptr<Transformer> trans_spmt_calib = std::make_shared<CalibrationTransformer>(
        "CalibrationTransformer", RecPmtType::PMT_3INCH, 13.0
    );
    m_pipe->addStep(trans_spmt_calib);

    std::shared_ptr<Transformer> trans_hama_q = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_20INCH_HAMAMATSU, 5.0
    );
    m_pipe->addStep(trans_hama_q);

    std::shared_ptr<Transformer> trans_nnvt_q = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_20INCH_NNVT, 20.0
    );
    m_pipe->addStep(trans_nnvt_q);
    
    std::shared_ptr<Transformer> trans = std::make_shared<WaterPhaseTransformer>(
        "WaterPhaseTransformer", 200, 0.0, 1000.0, 30.0, 0.2, 1500.0, 10u, 35.0
    );
    m_pipe->addStep(trans);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        std::make_shared<CdFhtPredictor<ParamsType::SingleCd>>(
            std::make_shared<WaterPhaseCdFht<ParamsType::SingleCd>>()
        ),
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti = std::make_shared<FhtMinimizerEstimator>("FhtMinimizerEstimator", g_opti, cost);
    m_pipe->addStep(esti);
}

void CdWaterPhaseStrategy::prepare() {
    getDefaultParams<ParamsType::SingleCd>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdWaterPhaseStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<ParamsType::SingleCd>::set(params, t_start, start, dir, length);
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
// =================================== CD WP Water Phase Strategy =================================
// ################################################################################################

void CdWpWaterPhaseStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<WaterPhaseInitializer>(
        "WaterPhaseInitializer", 200.0, 20.0, 0.25
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);
    m_pipe->addStep(g_trans_calib_wp);

    std::shared_ptr<Transformer> trans_hama_calib = std::make_shared<CalibrationTransformer>(
        "CalibrationTransformer", RecPmtType::PMT_20INCH_HAMAMATSU, -7.0
    );
    m_pipe->addStep(trans_hama_calib);
    
    std::shared_ptr<Transformer> trans_spmt_calib = std::make_shared<CalibrationTransformer>(
        "CalibrationTransformer", RecPmtType::PMT_3INCH, 13.0
    );
    m_pipe->addStep(trans_spmt_calib);

    std::shared_ptr<Transformer> trans_hama_q = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_20INCH_HAMAMATSU, 5.0
    );
    m_pipe->addStep(trans_hama_q);

    std::shared_ptr<Transformer> trans_nnvt_q = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_20INCH_NNVT, 20.0
    );
    m_pipe->addStep(trans_nnvt_q);
    
    std::shared_ptr<Transformer> trans = std::make_shared<WaterPhaseTransformer>(
        "WaterPhaseTransformer", 200, 0.0, 1000.0, 30.0, 0.2, 1500.0, 10u, 35.0
    );
    m_pipe->addStep(trans);

    std::shared_ptr<Transformer> trans_wp_clb = std::make_shared<CalibrationTransformer>(
        "CalibrationTransformer", RecPmtType::PMT_WP, 20.0
    );
    m_pipe->addStep(trans_wp_clb);

    std::shared_ptr<Transformer> trans_wp_q = std::make_shared<FhtChargeTholdTransformer>(
        "FhtChargeTholdTransformer", RecPmtType::PMT_WP, 5.0
    );
    m_pipe->addStep(trans_wp_q);

    std::shared_ptr<Transformer> trans_wp_geomtime = std::make_shared<WpGeomTimeTransformer>(
        "WpGeomTimeTransformer", 50.0, 100.0, 0.4, 0.8, 10.0, 0.05
    );
    m_pipe->addStep(trans_wp_geomtime);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        std::make_shared<CdWpFhtPredictor<ParamsType::SingleCd>>(
            std::make_shared<WaterPhaseCdFht<ParamsType::SingleCd>>(),
            std::make_shared<NoGeomWpFht<ParamsType::SingleCd>>()
        ),
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti = std::make_shared<FhtMinimizerEstimator>("FhtMinimizerEstimator", g_opti, cost);
    m_pipe->addStep(esti);
}

void CdWpWaterPhaseStrategy::prepare() {
    getDefaultParams<ParamsType::SingleCd>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdWpWaterPhaseStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<ParamsType::SingleCd>::set(params, t_start, start, dir, length);
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