#include "CdWpTtChi2RecTool.hpp"

#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(CdWpTtChi2RecTool);

CdWpTtChi2RecTool::CdWpTtChi2RecTool(const std::string& name) :
    ToolBase(name),
    m_fact{"StractegyFactory"},
    m_strat{nullptr},
    m_timer{name}
{
    declProp("WaterPhase", m_waterphase);
}

bool CdWpTtChi2RecTool::initialize() {
    if (m_waterphase) {
        m_fact.addStrategy(std::make_shared<CdWaterPhaseStrategy>());
        m_fact.addStrategy(std::make_shared<CdWpWaterPhaseStrategy>());
    }
    else {
        m_fact.addStrategy(std::make_shared<CdStrategy>());
        // m_fact.addStrategy(std::make_shared<CdStoppingStrategy>());
        // m_fact.addStrategy(std::make_shared<CdDoubleStrategy>());
        m_fact.addStrategy(std::make_shared<TtStrategy>());
        m_fact.addStrategy(std::make_shared<CdWpStrategy>());
        m_fact.addStrategy(std::make_shared<CdTtStrategy>());
        m_fact.addStrategy(std::make_shared<CdWpTtStrategy>());
    }
    if (!m_fact.initialize()) return false;
    LogInfo << m_name << " initialized successfully\n";
}

bool CdWpTtChi2RecTool::configure(const Params* params, const PmtTable* table) {
    m_params = params;
    m_ref_table = table;
    if (!initialize()) return false;
    LogInfo << m_name << " configured successfully\n";
    return true;
}

bool CdWpTtChi2RecTool::reconstruct(RecTrks* trks) {
    timer_guard tg(m_timer); // start the timer, and stop it when `tg` goes out of scope and if timer::stop() is not called

    TableConverter::convert(m_ref_table, m_table);
    double tot_pe = TableConverter::getTotPE();

    LogInfo << "Start reconstruction entry with " << m_table.size() << " PMTs, and " << tot_pe << " PEs\n";

    m_start = m_fact.construct(m_table);
    if (!m_strat) {
        LogInfo << "Skip this event due to no strategy found for this type of event\n";
        return true;
    }

    if (!m_strat->get()->estimate(m_table)) {
        LogInfo << "Skip this event due to failed estimation\n";
        return true;
    }

    std::vector<double> res = m_strat->get()->getFVars();
    double chi2 = m_strat->get()->getCost();
    m_strat->saveTrack(trks, res, chi2, tot_pe);
    
    return true;
}

bool CdWpTtChi2RecTool::finalize() {
    if (!m_fact.finalize()) return false;
    LogInfo << m_name << " finalized successfully\n";
    return true;
}