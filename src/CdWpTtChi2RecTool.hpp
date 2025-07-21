#ifndef CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_
#define CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_

#include "SniperKernel/ToolBase.h"
#include "RecTools/IRecMuonTool.h"

#include <any>
#include <memory>
#include <unordered_map>

#include "strategy/Factory.hpp"
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

    const Params* m_params;
    const PmtTable* m_ref_table;

    RecPmtTable m_table;

    bool m_auto_fact;
    bool m_water_phase;
    int m_manual_reco_mode;

    StrategyRegistry m_reg;
    std::string m_config_file;
    std::shared_ptr<StrategyFactory> m_fact;

    timer m_timer;

};

#endif // CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_