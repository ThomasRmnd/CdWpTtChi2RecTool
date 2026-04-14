#include "strategy/Strategy.hpp"

#include "chi2/fht/FhtChi2.hpp"
#include "chi2/tt/TtChi2.hpp"

#include "cost/fht/FhtCostFunction.hpp"
#include "cost/fht_tt/FhtTtCostFunction.hpp"
#include "cost/tt/TtCostFunction.hpp"

#include "estimator/ParamsModifierEstimator.hpp"
#include "estimator/fht/CorrectionMapLoopEstimator.hpp"
#include "estimator/fht/FhtMinimizerEstimator.hpp"
#include "estimator/fht/map/CorrectionMap3d.hpp"
#include "estimator/fht/map/CorrParam.hpp"
#include "estimator/fht/map/WpTimeShiftCorrectionMap.hpp"
#include "estimator/fht_tt/FhtTtCorrMapEstimator.hpp"
#include "estimator/optimizer/RootOptimizer.hpp"
#include "estimator/tt/TtMinimizerEstimator.hpp"
#include "estimator/tt/combinator/MaskHeightCartesianProdCombinator.hpp"
#include "estimator/tt/converter/FuzeNeighborConverter.hpp"

// #include "initializer/fht/ClusterBundleInitializer.hpp" WIP
#include "initializer/fht/ClusterMaxChargeInitializer.hpp"
#include "initializer/fht/WaterPhaseInitializer.hpp"
#include "initializer/tt/MinMaxHeightInitializer.hpp"

#include "predictor/fht/CdFhtPredictor.hpp"
#include "predictor/fht/CdWpFhtPredictor.hpp"
#include "predictor/fht/WpFhtPredictor.hpp"
#include "predictor/fht/fht/CdFht.hpp"
#include "predictor/fht/fht/WpFht.hpp"
#include "predictor/tt/TtPredictor.hpp"

#include "transformer/CalibTimeTransformer.hpp"
#include "transformer/EarlyLateFhtTransformer.hpp"
#include "transformer/FhtChargeTholdTransformer.hpp"
#include "transformer/PmtTypeTransformer.hpp"
#include "transformer/TtCrossTalkTransformer.hpp"
#include "transformer/WaterPhaseTransformer.hpp"

#define DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(Base, name, Derived, ...) \
    std::shared_ptr<Base> g_##name##_single = std::make_shared<Derived<SingleAcrylicParamsTag>>(__VA_ARGS__); \
    std::shared_ptr<Base> g_##name##_single_stopping = std::make_shared<Derived<SingleStoppingAcrylicParamsTag>>(__VA_ARGS__); \
    std::shared_ptr<Base> g_##name##_double = std::make_shared<Derived<DoubleAcrylicParamsTag>>(__VA_ARGS__);

#define DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(Base, name, Derived, ...) \
    std::shared_ptr<Base<SingleAcrylicParamsTag>> g_##name##_single = std::make_shared<Derived<SingleAcrylicParamsTag>>(__VA_ARGS__); \
    std::shared_ptr<Base<SingleStoppingAcrylicParamsTag>> g_##name##_single_stopping = std::make_shared<Derived<SingleStoppingAcrylicParamsTag>>(__VA_ARGS__); \
    std::shared_ptr<Base<DoubleAcrylicParamsTag>> g_##name##_double = std::make_shared<Derived<DoubleAcrylicParamsTag>>(__VA_ARGS__);

#define DEFINIT_GLOBAL_PREDICTOR(Base, name, Derived, fhtname) \
    std::shared_ptr<Base> g_##name##_single = std::make_shared<Derived<SingleAcrylicParamsTag>>(g_##fhtname##_single); \
    std::shared_ptr<Base> g_##name##_single_stopping = std::make_shared<Derived<SingleStoppingAcrylicParamsTag>>(g_##fhtname##_single_stopping); \
    std::shared_ptr<Base> g_##name##_double = std::make_shared<Derived<DoubleAcrylicParamsTag>>(g_##fhtname##_double);

#define DEFINIT_GLOBAL_PREDICTOR_CDWP(Base, name, Derived, cdfhtname, wpfhtname) \
    std::shared_ptr<Base> g_##name##_single = std::make_shared<Derived<SingleAcrylicParamsTag>>(g_##cdfhtname##_single, g_##wpfhtname##_single); \
    std::shared_ptr<Base> g_##name##_single_stopping = std::make_shared<Derived<SingleStoppingAcrylicParamsTag>>(g_##cdfhtname##_single_stopping, g_##wpfhtname##_single_stopping); \
    std::shared_ptr<Base> g_##name##_double = std::make_shared<Derived<DoubleAcrylicParamsTag>>(g_##cdfhtname##_double, g_##wpfhtname##_double);

std::shared_ptr<Chi2<FhtMethodTag>> g_chi2_fht = std::make_shared<FhtChi2>("FhtChi2");
std::shared_ptr<Chi2<TtMethodTag>> g_chi2_tt = std::make_shared<TtChi2>("TtChi2", 13.0);
std::shared_ptr<Chi2<TtMethodTag>> g_chi2_tt_joint = std::make_shared<TtChi2>("TtChi2_Joint", 130.0);

std::shared_ptr<CorrParam> g_corr_param_dist_proj_pmt_to_orig = std::make_shared<DistProjPmtToOrigCorrParam>();
std::shared_ptr<CorrParam> g_corr_param_dist_proj_pmt = std::make_shared<DistProjPmtCorrParam>();
std::shared_ptr<CorrParam> g_corr_param_angle = std::make_shared<AngleCorrParam>();
std::shared_ptr<CorrParam> g_corr_param_dist_track_to_center_squared = std::make_shared<DistTrackToCenterSquaredCorrParam>();

DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(CorrectionMap, corr_map_nnvt, CorrectionMap3d, 
    "CorrectionMap3d_NNVT", RecPmtType::PMT_20INCH_NNVT | RecPmtType::PMT_20INCH_HIGHQENNVT,
    "data/Reconstruction/RecMuon/CdWpTtChi2RecTool/RUN.9789-10084.reprod25c.correctionmap.root", "NNVTPmtsCorrectionMap_FhtChargeThold20",
    g_corr_param_dist_proj_pmt, g_corr_param_angle, g_corr_param_dist_track_to_center_squared
)
DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(CorrectionMap, corr_map_hamamatsu, CorrectionMap3d,
    "CorrectionMap3d_Hamamatsu", RecPmtType::PMT_20INCH_HAMAMATSU,
    "data/Reconstruction/RecMuon/CdWpTtChi2RecTool/RUN.9789-10084.reprod25c.correctionmap.root", "HamamatsuPmtsCorrectionMap_FhtChargeThold20",
    g_corr_param_dist_proj_pmt, g_corr_param_angle, g_corr_param_dist_track_to_center_squared
)
DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(CorrectionMap, corr_map_3inch, CorrectionMap3d,
    "CorrectionMap3d_3inch", RecPmtType::PMT_3INCH,
    "data/Reconstruction/RecMuon/CdWpTtChi2RecTool/RUN.9789-10084.reprod25c.correctionmap.root", "3inchPmtsCorrectionMap_EarlyLateFht135",
    g_corr_param_dist_proj_pmt, g_corr_param_angle, g_corr_param_dist_track_to_center_squared
)
DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(CorrectionMap, corr_map_wp, WpTimeShiftCorrectionMap,
    "WpTimeShiftCorrectionMap",
    "data/Reconstruction/RecMuon/CdWpTtChi2RecTool/RUN.9789-10084.reprod25c.correctionmap.root", "WpPmtsCorrectionMap_TimeGeom"
)

std::shared_ptr<Optimizer> g_opti = std::make_shared<RootOptimizer>(1000000, 100000, 0.001);

std::shared_ptr<Combinator> g_tt_combinator = std::make_shared<MaskHeightCartesianProdCombinator>("MaskHeightCartesianProdCombinator", 6);

std::shared_ptr<Converter> g_tt_converter = std::make_shared<FuzeNeighborConverter>("FuzeNeighborConverter");

std::shared_ptr<Initializer<TtMethodTag>> g_tt_initializer_minmax_height = std::make_shared<MinMaxHeightInitializer>("MinMaxHeightInitializer");

DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CdFht, cd_fht, NoRefractionLsCdFht)
DEFINIT_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(WpFht, wp_fht, NoDiffusionWpFht)

DEFINIT_GLOBAL_PREDICTOR(Predictor<FhtMethodTag>, pred_fht_cd_no_refr_ls, CdFhtPredictor, cd_fht)
DEFINIT_GLOBAL_PREDICTOR(Predictor<FhtMethodTag>, pred_fht_wp_no_hit, WpFhtPredictor, wp_fht)
DEFINIT_GLOBAL_PREDICTOR_CDWP(Predictor<FhtMethodTag>, pred_fht_no_refr_ls_no_hit, CdWpFhtPredictor, cd_fht, wp_fht)
DEFINIT_GLOBAL_BASED_ON_TRACK_PARAMS(Predictor<TtMethodTag>, pred_tt, TtPredictor)
std::shared_ptr<Predictor<TtMethodTag>> g_pred_tt_single_tt = std::make_shared<TtPredictor<SingleTtParamsTag>>();

std::shared_ptr<Transformer> g_trans_calib_hama = std::make_shared<CalibTimeTransformer>("CalibTimeTransformer_Hamamatsu", RecPmtType::PMT_20INCH_HAMAMATSU, 0.0);
std::shared_ptr<Transformer> g_trans_calib_nnvt = std::make_shared<CalibTimeTransformer>("CalibTimeTransformer_NNVT", RecPmtType::PMT_20INCH_NNVT, 0.0);
std::shared_ptr<Transformer> g_trans_calib_highqe = std::make_shared<CalibTimeTransformer>("CalibTimeTransformer_HighQENNVT", RecPmtType::PMT_20INCH_HIGHQENNVT, 0.0);
std::shared_ptr<Transformer> g_trans_calib_spmt = std::make_shared<CalibTimeTransformer>("CalibTimeTransformer_3inch", RecPmtType::PMT_3INCH, 0.0);
std::shared_ptr<Transformer> g_trans_calib_wp = std::make_shared<CalibTimeTransformer>("CalibTimeTransformer_WP", RecPmtType::PMT_WP, 0.0);
std::shared_ptr<Transformer> g_trans_tt_cross_talk = std::make_shared<TtCrossTalkTransformer>("TtCrossTalkTransformer");

// ################################################################################################
// ========================================= CD Strategy ==========================================
// ################################################################################################

void CdStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("CdStrategy__Pipeline");

    // ======================================= Initializer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "CdStrategy__ClusterMaxChargeInitializer", 1000, 200, 0.0, 1000.0, -10.0, 10.0, 20.0, 0.0, 50.0, 0.5, ClusterMaxChargeInitializer::Mode::CD_ONLY
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "CdStrategy__FhtChargeTholdTransformer_20inch", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    // std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
    //     "CdStrategy__EarlyLateFhtTransformer_3inch", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    // );
    std::shared_ptr<Transformer> trans_type_spmt = std::make_shared<PmtTypeTransformer>(
        "CdStrategy__PmtTypeTransformer_3inch", RecPmtType::PMT_3INCH
    );
    m_pipe->addStep(trans_type_spmt);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        "CdStrategy__FhtCostFunction",
        g_pred_fht_cd_no_refr_ls_single,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "CdStrategy__FhtMinimizerEstimator", g_opti, cost
    );
    // m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<CorrectionMapLoopEstimator>(
        "CdStrategy__CorrectionMapLoopEstimator", esti_min, 
        std::vector<std::shared_ptr<CorrectionMap>>{g_corr_map_nnvt_single, g_corr_map_hamamatsu_single, g_corr_map_3inch_single}, 2ul
    );
    // m_pipe->addStep(esti_corrmap_loop);
}

void CdStrategy::prepare() {
    getDefaultParams<SingleAcrylicParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<SingleAcrylicParamsTag>::set(params, t_start, start, dir, length);
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
    m_pipe = std::make_shared<Pipeline>("CdStoppingStrategy__Pipeline");

    // ======================================= Initializer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "CdStoppingStrategy__ClusterMaxChargeInitializer", 1000, 200, 0.0, 1000.0, -10.0, 10.0, 20.0, 0.0, 50.0, 0.5, ClusterMaxChargeInitializer::Mode::CD_ONLY
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "CdStoppingStrategy__FhtChargeTholdTransformer_20inch", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    // std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
    //     "CdStoppingStrategy__EarlyLateFhtTransformer_3inch", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    // );
    std::shared_ptr<Transformer> trans_type_spmt = std::make_shared<PmtTypeTransformer>(
        "CdStoppingStrategy__PmtTypeTransformer_3inch", RecPmtType::PMT_3INCH
    );
    m_pipe->addStep(trans_type_spmt);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        "CdStoppingStrategy__FhtCostFunction",
        g_pred_fht_cd_no_refr_ls_single,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "CdStoppingStrategy__FhtMinimizerEstimator", g_opti, cost
    );
    m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<CorrectionMapLoopEstimator>(
        "CdStoppingStrategy__CorrectionMapLoopEstimator", esti_min, 
        std::vector<std::shared_ptr<CorrectionMap>>{g_corr_map_nnvt_single, g_corr_map_hamamatsu_single, g_corr_map_3inch_single}, 2ul
    );
    m_pipe->addStep(esti_corrmap_loop);

}

void CdStoppingStrategy::prepare() {
    getDefaultParams<SingleAcrylicParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdStoppingStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<SingleAcrylicParamsTag>::set(params, t_start, start, dir, length);
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

/* WIP

void CdDoubleStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("CdDoubleStrategy__Pipeline");

    // ======================================= Initializer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterBundleInitializer>(
        "CdDoubleStrategy__ClusterBundleInitializer", 1000.0, 5000.0
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "CdDoubleStrategy__FhtChargeTholdTransformer_20inch", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    // std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
    //     "CdDoubleStrategy__EarlyLateFhtTransformer_3inch", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    // );
    std::shared_ptr<Transformer> trans_type_spmt = std::make_shared<PmtTypeTransformer>(
        "CdDoubleStrategy__PmtTypeTransformer_3inch", RecPmtType::PMT_3INCH
    );
    m_pipe->addStep(trans_type_spmt);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        "CdDoubleStrategy__FhtCostFunction",
        g_pred_fht_cd_no_refr_ls_double,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "CdDoubleStrategy__FhtMinimizerEstimator", g_opti, cost
    );
    m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<ParamsConstrainerEstimator<DoubleAcrylicParamsTag>>(
        "CdDoubleStrategy__ParamsConstrainerEstimator",
        std::make_shared<CorrectionMapLoopEstimator>(
            "CdDoubleStrategy__CorrectionMapLoopEstimator", esti_min, 
            std::vector<std::shared_ptr<CorrectionMap>>{g_corr_map_nnvt_double, g_corr_map_hamamatsu_double, g_corr_map_3inch_double}, 2ul
        ),
        std::array<bool, 8>{false, false, false, false, false, false, true, true}
    );
    m_pipe->addStep(esti_corrmap_loop);
}

void CdDoubleStrategy::prepare() {
    getDefaultParams<DoubleAcrylicParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdDoubleStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start_1, start_2, dir_1, dir_2, end_1, end_2;
    double t_start_1, t_start_2, length_1, length_2, t_end_1, t_end_2;
    TrackSetterHelper<DoubleAcrylicParamsTag>::set(params, t_start_1, start_1, dir_1, length_1);
    TrackSetterHelper<DoubleAcrylicParamsTag>::set(params, t_start_2, start_2, dir_2, length_2);
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

*/

// ################################################################################################
// ========================================= TT Strategy ==========================================
// ################################################################################################

void TtStrategy::create() {
    m_pipe = std::make_shared<Pipeline>("TtStrategy__Pipeline");

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_tt_cross_talk);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<Estimator> esti_tt = std::make_shared<TtMinimizerEstimator>(
        "TtStrategy__TtMinimizerEstimator",
        g_opti,
        std::make_shared<TtCostFunction>(
            "TtStrategy__TtCostFunction",
            g_pred_tt_single_tt,
            g_chi2_tt
        ),
        g_tt_converter, 25, g_tt_combinator, g_tt_initializer_minmax_height
    );
    m_pipe->addStep(esti_tt);
}

void TtStrategy::prepare() {
    getDefaultParams<SingleTtParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void TtStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<SingleTtParamsTag>::set(params, t_start, start, dir, length);
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
    m_pipe = std::make_shared<Pipeline>("CdWpStrategy__Pipeline");

    // ======================================= Initializer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);
    m_pipe->addStep(g_trans_calib_wp);

    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "CdWpStrategy__ClusterMaxChargeInitializer", 1000, 200, 0.0, 1000.0, -10.0, 10.0, 20.0, 0.0, 50.0, 0.5, ClusterMaxChargeInitializer::Mode::CDWP_COMBINED
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "CdWpStrategy__FhtChargeTholdTransformer_20inch", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    // std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
    //     "CdWpStrategy__EarlyLateFhtTransformer_3inch", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    // );
    std::shared_ptr<Transformer> trans_type_spmt = std::make_shared<PmtTypeTransformer>(
        "CdWpStrategy__PmtTypeTransformer_3inch", RecPmtType::PMT_3INCH
    );
    m_pipe->addStep(trans_type_spmt);

    std::shared_ptr<Transformer> trans_q_wp = std::make_shared<FhtChargeTholdTransformer>(
        "CdWpStrategy__FhtChargeTholdTransformer_WP", RecPmtType::PMT_WP, 30.0
    );
    m_pipe->addStep(trans_q_wp);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        "CdWpStrategy__FhtCostFunction",
        g_pred_fht_no_refr_ls_no_hit_single,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "CdWpStrategy__FhtMinimizerEstimator", g_opti, cost
    );
    // m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<CorrectionMapLoopEstimator>(
        "CdWpStrategy__CorrectionMapLoopEstimator", esti_min, 
        std::vector<std::shared_ptr<CorrectionMap>>{g_corr_map_nnvt_single, g_corr_map_hamamatsu_single, g_corr_map_3inch_single}, 2ul
    );
    // m_pipe->addStep(esti_corrmap_loop);
}

void CdWpStrategy::prepare() {
    getDefaultParams<SingleAcrylicParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdWpStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<SingleAcrylicParamsTag>::set(params, t_start, start, dir, length);
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
    m_pipe = std::make_shared<Pipeline>("CdTtStrategy__Pipeline");

    // ======================================= Initializer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "CdTtStrategy__ClusterMaxChargeInitializer", 1000, 200, 0.0, 1000.0, -10.0, 10.0, 20.0, 0.0, 50.0, 0.5, ClusterMaxChargeInitializer::Mode::CD_ONLY
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_tt_cross_talk);

    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "CdTtStrategy__FhtChargeTholdTransformer", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    // std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
    //     "CdTtStrategy__EarlyLateFhtTransformer_3inch", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    // );
    std::shared_ptr<Transformer> trans_type_spmt = std::make_shared<PmtTypeTransformer>(
        "CdTtStrategy__PmtTypeTransformer_3inch", RecPmtType::PMT_3INCH
    );
    m_pipe->addStep(trans_type_spmt);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        "CdTtStrategy__FhtMethodTag",
        g_pred_fht_cd_no_refr_ls_single,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "CdTtStrategy__FhtMinimizerEstimator", g_opti, cost
    );
    m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<CorrectionMapLoopEstimator>(
        "CdTtStrategy__CorrectionMapLoopEstimator", esti_min, 
        std::vector<std::shared_ptr<CorrectionMap>>{g_corr_map_nnvt_single, g_corr_map_hamamatsu_single, g_corr_map_3inch_single}, 2ul
    );

    std::shared_ptr<CostFunction<TtMethodTag>> tt_cost = std::make_shared<TtCostFunction>(
        "CdTtStrategy__TtCostFunction",
        std::make_shared<TtPredictor<SingleTtParamsTag>>(),
        g_chi2_tt
    );

    std::shared_ptr<Estimator> esti_tt = std::make_shared<TtMinimizerEstimator>(
        "CdTtStrategy__TtMinimizerEstimator",
        g_opti, tt_cost,
        g_tt_converter, 25, g_tt_combinator, g_tt_initializer_minmax_height
    );

    std::shared_ptr<CostFunction<FhtTtMethodTag>> fht_tt_cost = std::make_shared<FhtTtCostFunction>(
        "CdTtStrategy__FhtTtCostFunction",
        g_pred_fht_cd_no_refr_ls_single,
        g_chi2_fht,
        g_pred_tt_single,
        g_chi2_tt_joint
    );

    std::shared_ptr<Estimator> esti_fhttt_corrmap = std::make_shared<FhtTtCorrMapEstimator>(
        "CdTtStrategy__FhtTtCorrectionMapEstimator",
        g_opti, fht_tt_cost, 
        std::dynamic_pointer_cast<CorrectionMapLoopEstimator>(esti_corrmap_loop), 
        std::dynamic_pointer_cast<TtMinimizerEstimator>(esti_tt),
        10, 2
    );
    m_pipe->addStep(esti_fhttt_corrmap);
}

void CdTtStrategy::prepare() {
    getDefaultParams<SingleAcrylicParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdTtStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<SingleAcrylicParamsTag>::set(params, t_start, start, dir, length);
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
    m_pipe = std::make_shared<Pipeline>("CdWpTtStrategy__Pipeline");

    // ======================================= Initializer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);
    m_pipe->addStep(g_trans_calib_wp);

    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<ClusterMaxChargeInitializer>(
        "CdWpTtStrategy__ClusterMaxChargeInitializer", 1000, 200, 0.0, 1000.0, -10.0, 10.0, 20.0, 0.0, 50.0, 0.5, ClusterMaxChargeInitializer::Mode::CDWP_COMBINED
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_tt_cross_talk);

    std::shared_ptr<Transformer> trans_q_lpmt = std::make_shared<FhtChargeTholdTransformer>(
        "CdWpTtStrategy__FhtChargeTholdTransformer", RecPmtType::PMT_20INCH, 20.0
    );
    m_pipe->addStep(trans_q_lpmt);
    
    // std::shared_ptr<Transformer> trans_fht_spmt = std::make_shared<EarlyLateFhtTransformer>(
    //     "CdWpTtStrategy__EarlyLateFhtTransformer_3inch", RecPmtType::PMT_3INCH, 1000, 0.0, 1000.0, 2.0, 2.0, 135.0
    // );
    std::shared_ptr<Transformer> trans_type_spmt = std::make_shared<PmtTypeTransformer>(
        "CdWpTtStrategy__PmtTypeTransformer_3inch", RecPmtType::PMT_3INCH
    );
    m_pipe->addStep(trans_type_spmt);

    std::shared_ptr<Transformer> trans_q_wp = std::make_shared<FhtChargeTholdTransformer>(
        "CdWpTtStrategy__FhtChargeTholdTransformer", RecPmtType::PMT_WP, 30.0
    );
    m_pipe->addStep(trans_q_wp);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        "CdWpTtStrategy__FhtCostFunction",
        g_pred_fht_no_refr_ls_no_hit_single,
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti_min = std::make_shared<FhtMinimizerEstimator>(
        "CdWpTtStrategy__FhtMinimizerEstimator", g_opti, cost
    );
    m_pipe->addStep(esti_min);

    // ===================================== 2+ Minimization ======================================
    std::shared_ptr<Estimator> esti_corrmap_loop = std::make_shared<CorrectionMapLoopEstimator>(
        "CdWpTtStrategy__CorrectionMapLoopEstimator", esti_min, 
        std::vector<std::shared_ptr<CorrectionMap>>{g_corr_map_nnvt_single, g_corr_map_hamamatsu_single, g_corr_map_3inch_single}, 2ul
    );

    std::shared_ptr<CostFunction<TtMethodTag>> tt_cost = std::make_shared<TtCostFunction>(
        "CdWpTtStrategy__TtCostFunction",
        g_pred_tt_single_tt,
        g_chi2_tt
    );

    std::shared_ptr<Estimator> esti_tt = std::make_shared<TtMinimizerEstimator>(
        "CdWpTtStrategy__TtMinimizerEstimator",
        g_opti, tt_cost,
        g_tt_converter, 25, g_tt_combinator, g_tt_initializer_minmax_height
    );

    std::shared_ptr<CostFunction<FhtTtMethodTag>> fht_tt_cost = std::make_shared<FhtTtCostFunction>(
        "CdWpTtStrategy__FhtTtCostFunction",
        g_pred_fht_no_refr_ls_no_hit_single,
        g_chi2_fht,
        g_pred_tt_single,
        g_chi2_tt_joint
    );

    std::shared_ptr<Estimator> esti_fhttt_corrmap = std::make_shared<FhtTtCorrMapEstimator>(
        "CdWpTtStrategy__FhtTtCorrectionMapEstimator",
        g_opti, fht_tt_cost, 
        std::dynamic_pointer_cast<CorrectionMapLoopEstimator>(esti_corrmap_loop), 
        std::dynamic_pointer_cast<TtMinimizerEstimator>(esti_tt),
        10, 2
    );
    m_pipe->addStep(esti_fhttt_corrmap);
}

void CdWpTtStrategy::prepare() {
    getDefaultParams<SingleAcrylicParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdWpTtStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<SingleAcrylicParamsTag>::set(params, t_start, start, dir, length);
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
    m_pipe = std::make_shared<Pipeline>("CdWaterPhaseStrategy__Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<WaterPhaseInitializer>(
        "CdWaterPhaseStrategy__WaterPhaseInitializer", 200.0, 20.0, 0.25
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);

    std::shared_ptr<Transformer> trans_hama_q = std::make_shared<FhtChargeTholdTransformer>(
        "CdWaterPhaseStrategy__FhtChargeTholdTransformer", RecPmtType::PMT_20INCH_HAMAMATSU, 5.0
    );
    m_pipe->addStep(trans_hama_q);

    std::shared_ptr<Transformer> trans_nnvt_q = std::make_shared<FhtChargeTholdTransformer>(
        "CdWaterPhaseStrategy__FhtChargeTholdTransformer", RecPmtType::PMT_20INCH_NNVT, 20.0
    );
    m_pipe->addStep(trans_nnvt_q);
    
    std::shared_ptr<Transformer> trans = std::make_shared<WaterPhaseTransformer>(
        "CdWaterPhaseStrategy__WaterPhaseTransformer", 200, 0.0, 1000.0, 30.0, 0.2, 1500.0, 10u, 35.0
    );
    m_pipe->addStep(trans);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        "CdWaterPhaseStrategy__FhtCostFunction",
        std::make_shared<CdFhtPredictor<SingleCdParamsTag>>(
            std::make_shared<WaterPhaseCdFht<SingleCdParamsTag>>()
        ),
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti = std::make_shared<FhtMinimizerEstimator>("CdWaterPhaseStrategy__FhtMinimizerEstimator", g_opti, cost);
    m_pipe->addStep(esti);
}

void CdWaterPhaseStrategy::prepare() {
    getDefaultParams<SingleCdParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdWaterPhaseStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<SingleCdParamsTag>::set(params, t_start, start, dir, length);
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
    m_pipe = std::make_shared<Pipeline>("CdWpWaterPhaseStrategy__Pipeline");

    // ======================================= Initializer ========================================
    std::shared_ptr<Initializer<FhtMethodTag>> init = std::make_shared<WaterPhaseInitializer>(
        "CdWpWaterPhaseStrategy__WaterPhaseInitializer", 200.0, 20.0, 0.25
    );
    m_pipe->addStep(init);

    // ======================================= Transformer ========================================
    m_pipe->addStep(g_trans_calib_hama);
    m_pipe->addStep(g_trans_calib_highqe);
    m_pipe->addStep(g_trans_calib_nnvt);
    m_pipe->addStep(g_trans_calib_spmt);
    m_pipe->addStep(g_trans_calib_wp);

    std::shared_ptr<Transformer> trans_hama_q = std::make_shared<FhtChargeTholdTransformer>(
        "CdWpWaterPhaseStrategy__FhtChargeTholdTransformer", RecPmtType::PMT_20INCH_HAMAMATSU, 5.0
    );
    m_pipe->addStep(trans_hama_q);

    std::shared_ptr<Transformer> trans_nnvt_q = std::make_shared<FhtChargeTholdTransformer>(
        "CdWpWaterPhaseStrategy__FhtChargeTholdTransformer", RecPmtType::PMT_20INCH_NNVT, 20.0
    );
    m_pipe->addStep(trans_nnvt_q);
    
    std::shared_ptr<Transformer> trans = std::make_shared<WaterPhaseTransformer>(
        "CdWpWaterPhaseStrategy__WaterPhaseTransformer", 200, 0.0, 1000.0, 30.0, 0.2, 1500.0, 10u, 35.0
    );
    m_pipe->addStep(trans);

    std::shared_ptr<Transformer> trans_wp_q = std::make_shared<FhtChargeTholdTransformer>(
        "CdWpWaterPhaseStrategy__FhtChargeTholdTransformer", RecPmtType::PMT_WP, 30.0
    );
    m_pipe->addStep(trans_wp_q);

    // ===================================== 1st Minimization =====================================
    std::shared_ptr<CostFunction<FhtMethodTag>> cost = std::make_shared<FhtCostFunction>(
        "CdWpWaterPhaseStrategy__FhtCostFunction",
        std::make_shared<CdWpFhtPredictor<SingleCdParamsTag>>(
            std::make_shared<WaterPhaseCdFht<SingleCdParamsTag>>(),
            std::make_shared<NoDiffusionWpFht<SingleCdParamsTag>>()
        ),
        g_chi2_fht
    );

    std::shared_ptr<Estimator> esti = std::make_shared<FhtMinimizerEstimator>("CdWpWaterPhaseStrategy__FhtMinimizerEstimator", g_opti, cost);
    m_pipe->addStep(esti);
}

void CdWpWaterPhaseStrategy::prepare() {
    getDefaultParams<SingleCdParamsTag>();
    m_pipe->setParams(m_params, m_steps, m_names);
}

void CdWpWaterPhaseStrategy::save(RecTrks* tracks, double totpe) {
    const double* params = m_pipe->getParams().data();
    double cost = m_pipe->getCost();
    vec3 start, dir, end;
    double t_start, length, t_end;
    TrackSetterHelper<SingleCdParamsTag>::set(params, t_start, start, dir, length);
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