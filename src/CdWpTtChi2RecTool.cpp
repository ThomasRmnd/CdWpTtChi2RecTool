#include "CdWpTtChi2RecTool.hpp"

#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(CdWpTtChi2RecTool);

CdWpTtChi2RecTool::CdWpTtChi2RecTool(const std::string& name) :
    ToolBase(name),
    m_reg{"StrategyRegistry"},
    m_timer{name}
{
    declProp("UseAutomaticFactory", m_auto_fact = false);
    declProp("WaterPhase", m_water_phase = false);
    declProp("ManualReconstructionMode", m_manual_reco_mode = 1); // 1 = single through-going, 2 = single stopping, 3 = double
    declProp("ConfigMap", m_config_map);
}

bool CdWpTtChi2RecTool::initialize() {
    if (m_water_phase) {
        m_reg.book(std::make_shared<CdWaterPhaseStrategy>());
        m_reg.book(std::make_shared<CdWpWaterPhaseStrategy>());
    }
    else {
        m_reg.book(std::make_shared<CdStrategy>());
        m_reg.book(std::make_shared<CdStoppingStrategy>());
        m_reg.book(std::make_shared<CdDoubleStrategy>());
        m_reg.book(std::make_shared<TtStrategy>());
        m_reg.book(std::make_shared<CdWpStrategy>());
        m_reg.book(std::make_shared<CdTtStrategy>());
        m_reg.book(std::make_shared<CdWpTtStrategy>());
    }
    if (m_auto_fact) m_fact = std::make_shared<AutomaticStrategyFactory>(m_reg, nullptr);
    else {
        if (m_manual_reco_mode == 1) m_fact = std::make_shared<ManualStrategyFactory>(m_reg, ParamsType::SingleAcrylic);
        else if (m_manual_reco_mode == 2) m_fact = std::make_shared<ManualStrategyFactory>(m_reg, ParamsType::SingleStoppingAcrylic);
        else if (m_manual_reco_mode == 3) m_fact = std::make_shared<ManualStrategyFactory>(m_reg, ParamsType::DoubleAcrylic);
        else if (m_manual_reco_mode == 4) m_fact = std::make_shared<ManualStrategyFactory>(m_reg, ParamsType::SingleCd);
        else {
            LogError << "Reconstruction mode " << m_manual_reco_mode << " is not handled for the moment. Abording\n";
            return false;
        }
    }
    if (!m_reg.initialize()) return false;
    LogInfo << m_name << " initialized successfully\n";

    LogInfo << "Configuration Map:\n";
    for (const auto& [key, value_any] : m_config_map) {
        LogInfo << "  " << key << " = ";

        if (value_any.type() == typeid(int)) {
            LogInfo << std::any_cast<int>(value_any);
        }
        else if (value_any.type() == typeid(double)) {
            LogInfo << std::any_cast<double>(value_any);
        }
        else if (value_any.type() == typeid(std::string)) {
            LogInfo << std::any_cast<std::string>(value_any);
        }
        else if (value_any.type() == typeid(bool)) {
            LogInfo << (std::any_cast<bool>(value_any) ? "true" : "false");
        }
        else if (value_any.type() == typeid(RecPmtType)) {
            LogInfo << static_cast<int>(std::any_cast<RecPmtType>(value_any)); 
        }
        else {
            LogInfo << "<unprintable type: " << value_any.type().name() << ">";
        }
        LogInfo << "\n";
    }

    return true;
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
    double totpe = TableConverter::getTotPE();

    LogInfo << "Start reconstruction entry with " << m_table.size() << " PMTs, and " << totpe << " PEs\n";

    std::shared_ptr<Strategy> strat = m_fact->construct(m_table);
    if (!strat) {
        LogInfo << "No strategy found for this type of event. Skipping\n";
        return true;
    }
    strat->prepare();
    std::shared_ptr<Pipeline> pipe = strat->pipeline();
    if (!pipe) {
        LogError << "Pipeline is not set for this given strategy. Abording\n";
        return false;
    }

    if (!pipe->estimate(m_table)) {
        LogInfo << "Skip this event due to failed estimation\n";
        return true;
    }

    strat->save(trks, totpe);

    return true;
}

bool CdWpTtChi2RecTool::finalize() {
    if (!m_reg.finalize()) return false;
    LogInfo << m_name << " finalized successfully\n";
    return true;
}