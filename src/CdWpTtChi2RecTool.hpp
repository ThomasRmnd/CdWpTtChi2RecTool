#ifndef CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_
#define CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_

#include "SniperKernel/ToolBase.h"
#include "RecTools/IRecMuonTool.h"

#include <memory>

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

    bool m_waterphase;
    StrategyFactory m_fact;
    std::shared_ptr<Strategy> m_strat;

    timer m_timer;

};

#endif // CDWPTTCHI2RECTOOL_CDWPTTCHI2RECTOOL_HPP_