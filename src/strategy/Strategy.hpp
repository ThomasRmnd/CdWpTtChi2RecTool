#ifndef CDWPTTCHI2RECTOOL_STRATEGY_STRATEGY_HPP_
#define CDWPTTCHI2RECTOOL_STRATEGY_STRATEGY_HPP_

#include "SniperKernel/ToolBase.h"

#include "RecTools/IRecMuonTool.h"

#include "estimator/Pipeline.hpp"
#include "utils/TrackParams.hpp"

enum class DetectorType {
    NONE = 0,
    CD = 1 << 0,
    WP = 1 << 1,
    TT = 1 << 2
};

inline DetectorType operator|(const DetectorType& a, const DetectorType& b) {
    return static_cast<DetectorType>(static_cast<int>(a) | static_cast<int>(b));
}

inline DetectorType operator&(const DetectorType& a, const DetectorType& b) {
    return static_cast<DetectorType>(static_cast<int>(a) & static_cast<int>(b));
}

inline DetectorType operator^(const DetectorType& a, const DetectorType& b) {
    return static_cast<DetectorType>(static_cast<int>(a) ^ static_cast<int>(b));
}

inline DetectorType operator~(const DetectorType& a) {
    return static_cast<DetectorType>(~static_cast<int>(a));
}

inline DetectorType& operator|=(DetectorType& a, const DetectorType& b) {
    return a = a | b;
}

inline DetectorType& operator&=(DetectorType& a, const DetectorType& b) {
    return a = a & b;
}

inline DetectorType& operator^=(DetectorType& a, const DetectorType& b) {
    return a = a ^ b;
}

typedef std::pair<ParamsType, DetectorType> StrategyType;

struct StrategyType {

    ParamsType params;
    DetectorType detector;

};

inline bool operator==(const StrategyType& lhs, const StrategyType& rhs) {
    return (lhs.params == rhs.params) && (lhs.detector == rhs.detector);
}

namespace std {
template<>
struct hash<StrategyType> {
    std::size_t operator()(const StrategyType& key) const {
        return (static_cast<std::size_t>(key.params) << 3) ^ static_cast<std::size_t>(key.detector); //  << 3, since 3 detectors
    }
};
} // namespace std

class Strategy {

public:

    const StrategyType type;

    Strategy(const StrategyType& type_) : type(type_) {};

    virtual ~Strategy() = default;

    std::shared_ptr<Pipeline> pipeline() const { return m_pipe; };

    virtual void create() = 0;
    virtual void setDefaultParams() = 0;
    virtual void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) = 0;

protected:

    std::shared_ptr<Pipeline> m_pipe;

    std::vector<double> m_default_params;
    std::vector<double> m_steps;
    std::vector<std::string> m_names;

    template<typename TParams>
    void getParams() {
        m_default_params = std::vector<double>(TParams::default_params, TParams::default_params + TParams::nb_params);
        m_steps = std::vector<double>(TParams::steps, TParams::steps + TParams::nb_params);
        m_names = std::vector<std::string>(TParams::names, TParams::names + TParams::nb_params);
    }

};

class CdStrategy : public Strategy {

public:

    CdStrategy() : Strategy({ParamsType::SingleAcrylic, DetectorType::CD}) {}
    ~CdStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdStoppingStrategy : public Strategy {

public:

    CdStoppingStrategy() : Strategy({ParamsType::SingleStoppingAcrylic, DetectorType::CD}) {}
    ~CdStoppingStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdDoubleStrategy : public Strategy {

public:

    CdDoubleStrategy() : Strategy({ParamsType::DoubleAcrylic, DetectorType::CD}) {}
    ~CdDoubleStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class TtStrategy : public Strategy {

public:

    TtStrategy() : Strategy({ParamsType::SingleTt, DetectorType::TT}) {}
    ~TtStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdWpStrategy : public Strategy {

public:

    CdWpStrategy() : Strategy({ParamsType::SingleAcrylic, DetectorType::CD | DetectorType::WP}) {}
    ~CdWpStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdTtStrategy : public Strategy {

public:

    CdTtStrategy() : Strategy({ParamsType::SingleAcrylic, DetectorType::CD | DetectorType::TT}) {}
    ~CdTtStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdWpTtStrategy : public Strategy {

public:

    CdWpTtStrategy() : Strategy({ParamsType::SingleAcrylic, DetectorType::CD | DetectorType::WP | DetectorType::TT}) {}
    ~CdWpTtStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdWaterPhaseStrategy : public Strategy {

public:

    CdWaterPhaseStrategy() : Strategy({ParamsType::SingleCd, DetectorType::CD}) {}
    ~CdWaterPhaseStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdWpWaterPhaseStrategy : public Strategy {

public:

    CdWpWaterPhaseStrategy() : Strategy({ParamsType::SingleCd, DetectorType::CD | DetectorType::WP}) {}
    ~CdWpWaterPhaseStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

#define DECLARE_GLOBAL_BASED_ON_TRACK_PARAMS(Class, name) \
    extern std::shared_ptr<Class> g_##name##_single; \
    extern std::shared_ptr<Class> g_##name##_single_stopping; \
    extern std::shared_ptr<Class> g_##name##_double;

#define DECLARE_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(Class, name) \
    extern std::shared_ptr<Class<ParamsType::SingleAcrylic>>         g_##name##_single; \
    extern std::shared_ptr<Class<ParamsType::SingleStoppingAcrylic>> g_##name##_single_stopping; \
    extern std::shared_ptr<Class<ParamsType::DoubleAcrylic>>         g_##name##_double;

#include <memory>

#include "chi2/Chi2.hpp"

extern std::shared_ptr<Chi2<FhtMethodTag>> g_chi2_fht;
extern std::shared_ptr<Chi2<TtMethodTag>> g_chi2_tt;
extern std::shared_ptr<Chi2<TtMethodTag>> g_chi2_tt_joint;

#include "estimator/fht/map/CorrParam.hpp"

extern std::shared_ptr<CorrParam> g_corr_param_dist_proj_pmt_to_orig;
extern std::shared_ptr<CorrParam> g_corr_param_angle;
extern std::shared_ptr<CorrParam> g_corr_param_dist_track_to_center_squared;

#include "estimator/fht/map/CorrectionMap.hpp"

DECLARE_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CorrectionMap, corr_map_nnvt)
DECLARE_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CorrectionMap, corr_map_hamamatsu)
DECLARE_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CorrectionMap, corr_map_3inch)
DECLARE_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CorrectionMap, corr_map_wp)

#include "estimator/optimizer/Optimizer.hpp"

extern std::shared_ptr<Optimizer> g_opti;

#include "predictor/Predictor.hpp"

DECLARE_GLOBAL_BASED_ON_TRACK_PARAMS(Predictor<FhtMethodTag>, pred_fht_cd_no_refl_ls)
DECLARE_GLOBAL_BASED_ON_TRACK_PARAMS(Predictor<FhtMethodTag>, pred_fht_wp_no_hit)
DECLARE_GLOBAL_BASED_ON_TRACK_PARAMS(Predictor<FhtMethodTag>, pred_fht_no_refl_ls_no_hit)
DECLARE_GLOBAL_BASED_ON_TRACK_PARAMS(Predictor<TtMethodTag>, pred_tt)

#include "predictor/fht/fht/CdFht.hpp"
#include "predictor/fht/fht/WpFht.hpp"

DECLARE_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(CdFht, cd_fht)
DECLARE_GLOBAL_BASED_ON_TEMPLATE_TRACK_PARAMS(WpFht, wp_fht)


#endif // CDWPTTCHI2RECTOOL_STRATEGY_STRATEGY_HPP_