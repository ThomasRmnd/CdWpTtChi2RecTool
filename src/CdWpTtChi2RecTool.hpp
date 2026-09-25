#ifndef CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_
#define CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_

#include <memory>

#include "Geometry/IPMTParamSvc.h"
#include "RecTools/IRecMuonTool.h"
#include "SniperKernel/ToolBase.h"

#include "strategy/Factory.hpp"
#include "trigger/TriggerTimeCorrelation.hpp"
#include "utils/RecPmtProp.hpp"
#include "utils/timer.hpp"

class CdWpTtChi2RecTool : public ToolBase, public IRecMuonTool {

public:

    CdWpTtChi2RecTool(const std::string& name);

    ~CdWpTtChi2RecTool() override = default;

    bool initialize() override;
    bool finalize() override;

    bool configure(const Params* params, const PmtTable* table) override;
    bool reconstruct(RecTrks* trks) override;

private:

    bool m_auto_fact;
    bool m_water_phase;
    int m_manual_reco_mode;
    std::string m_config_file;

    const Params* c_params;
    const PmtTable* c_ref_table;
    RecPmtTable m_table;

    StrategyRegistry m_reg;
    SniperJSON m_json;
    std::shared_ptr<StrategyFactory> m_fact;
    std::shared_ptr<TriggerTimeCorrelation> m_trig_corr;

    JM::NavBuffer* m_buf;
    IPMTParamSvc* m_pmt_svc;

    timer m_timer;

};

#endif // CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_