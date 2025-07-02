#ifndef CDWPTTCHI2RECTOOL_FACTORY_STRATEGY_HPP_
#define CDWPTTCHI2RECTOOL_FACTORY_STRATEGY_HPP_

#include "SniperKernel/ToolBase.h"

#include "RecTools/IRecMuonTool.h"

#include "estimator/Pipeline.hpp"

enum class StrategyType {
    NONE = 0,
    CD = 1 << 0,
    WP = 1 << 1,
    TT = 1 << 2
};

inline StrategyType operator|(const StrategyType& a, const StrategyType& b) {
    return static_cast<StrategyType>(static_cast<int>(a) | static_cast<int>(b));
}

inline StrategyType operator&(const StrategyType& a, const StrategyType& b) {
    return static_cast<StrategyType>(static_cast<int>(a) & static_cast<int>(b));
}

inline StrategyType operator^(const StrategyType& a, const StrategyType& b) {
    return static_cast<StrategyType>(static_cast<int>(a) ^ static_cast<int>(b));
}

inline StrategyType operator~(const StrategyType& a) {
    return static_cast<StrategyType>(~static_cast<int>(a));
}

inline StrategyType& operator|=(const StrategyType&& a, const StrategyType& b) {
    return a = a | b;
}

inline StrategyType& operator&=(const StrategyType&& a, const StrategyType& b) {
    return a = a & b;
}

inline StrategyType& operator^=(const StrategyType&& a, const StrategyType& b) {
    return a = a ^ b;
}

class Strategy {

public:

    StrategyType type;

    Strategy(const StrategyType& t) : type(t) {};

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

    CdStrategy() : Strategy(StrategyType::CD) {};
    ~CdStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdDoubleStrategy : public Strategy {

public:

    CdDoubleStrategy() : Strategy(StrategyType::CD) {};
    ~CdDoubleStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class TtStrategy : public Strategy {

public:

    TtStrategy() : Strategy(StrategyType::TT) {};
    ~TtStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdWpStrategy : public Strategy {

public:

    CdWpStrategy() : Strategy(StrategyType::CD | StrategyType::WP) {};
    ~CdWpStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdTtStrategy : public Strategy {

public:

    CdTtStrategy() : Strategy(StrategyType::CD | StrategyType::TT) {};
    ~CdTtStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdWpTtStrategy : public Strategy {

public:

    CdWpTtStrategy() : Strategy(StrategyType::CD | StrategyType::WP | StrategyType::TT) {};
    ~CdWpTtStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdWaterPhaseStrategy : public Strategy {

public:

    CdWaterPhaseStrategy() : Strategy(StrategyType::CD) {};
    ~CdWaterPhaseStrategy() override = default;

    void create() override;
    void setDefaultParams() override;
    void saveTrack(RecTrks* tracks, const std::vector<double>& params, double cost, double tot_pe) override;

};

class CdWpWaterPhaseStrategy : public Strategy {

public:

    CdWpWaterPhaseStrategy() : Strategy(StrategyType::CD | StrategyType::WP) {};
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

#include "estimator/fht/correction_map/CorrParam.hpp"

extern std::shared_ptr<CorrParam> g_corr_param_dist_proj_pmt_to_orig;
extern std::shared_ptr<CorrParam> g_corr_param_angle;
extern std::shared_ptr<CorrParam> g_corr_param_dist_track_to_center_squared;

#include "estimator/fht/correction_map/CorrectionMap.hpp"

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


#endif // CDWPTTCHI2RECTOOL_FACTORY_STRATEGY_HPP_