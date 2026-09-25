#include "CdWpTtChi2RecTool.hpp"

#include <fstream>

#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/SniperJSON.h"
#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"

#include "DataPathHelper/Path.hh"

DECLARE_TOOL(CdWpTtChi2RecTool);

CdWpTtChi2RecTool::CdWpTtChi2RecTool(const std::string& name) :
    ToolBase{name},
    m_reg{"StrategyRegistry"},
    m_timer{name}
{
    declProp("UseAutomaticFactory", m_auto_fact = false);
    declProp("WaterPhase", m_water_phase = false);
    declProp("ManualReconstructionMode", m_manual_reco_mode = 1); // 1 = single through-going, 2 = single stopping, 3 = double
    declProp("ConfigFile", m_config_file);
}

bool CdWpTtChi2RecTool::initialize() {
    if (m_water_phase) {
        m_reg.book(std::make_shared<CdWaterPhaseStrategy>());
        m_reg.book(std::make_shared<CdWpWaterPhaseStrategy>());
    }
    else {
        m_reg.book(std::make_shared<CdStrategy>());
        m_reg.book(std::make_shared<CdStoppingStrategy>());
        // m_reg.book(std::make_shared<CdDoubleStrategy>());
        m_reg.book(std::make_shared<TtStrategy>());
        m_reg.book(std::make_shared<CdWpStrategy>());
        m_reg.book(std::make_shared<CdTtStrategy>());
        m_reg.book(std::make_shared<CdWpTtStrategy>());
    }
    if (m_auto_fact) {
        std::shared_ptr<Classifier> clsfr = nullptr;
        m_fact = std::make_shared<AutomaticStrategyFactory>(m_reg, clsfr);
    }
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

    std::string filename = JUNO::Path::resolve(m_config_file);
    if (filename.empty()) {
        LogWarn << "Cannot resolve config file: " << m_config_file << ", using default configuration\n";
    }
    else {
        std::ifstream ifs(filename);
        m_json = SniperJSON::load(ifs);
    }

    m_reg.configure(m_json);
    if (!m_reg.initialize()) return false;

    m_trig_corr = std::make_shared<TriggerTimeCorrelation>("TriggerTimeCorrelation");
    if (!m_trig_corr) {
        LogError << "Cannot create trigger time correlation tool\n";
        return false;
    }

    SniperDataPtr<JM::NavBuffer> navbuf(getParent(), "/Event");
    if (navbuf.invalid()) {
        LogError << "Cannot get the NavBuffer @ /Event\n";
        return false;
    }
    m_buf = navbuf.data();

    SniperPtr<IPMTParamSvc> pmtsvc(*getRoot(), "PMTParamSvc");
    if (pmtsvc.invalid()) {
        LogError << "Cannot get the PMTParamSvc\n";
        return false;
    }
    m_pmt_svc = pmtsvc.data();
    TableConverter::setPMTSvc(m_pmt_svc);

    LogInfo << m_name << " initialized successfully\n";

    return true;
}

bool CdWpTtChi2RecTool::configure(const Params* params, const PmtTable* table) {
    c_params = params;
    c_ref_table = table;
    if (!initialize()) return false;
    LogInfo << m_name << " configured successfully\n";
    return true;
}

bool CdWpTtChi2RecTool::reconstruct(RecTrks* trks) {
    timer_guard tg(m_timer); // start the timer, and stop it when `tg` goes out of scope and if timer::stop() is not called

    

    TableConverter::convert(c_ref_table, m_table);
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

    const std::vector<double>& params = pipe->getParams();
    for (double par : params) {
        if (std::isnan(par)) {
            LogInfo << "One of the track parameters is NaN. Skipping\n";
            return true;
        }
    }

    strat->save(trks, totpe);

    return true;
}

bool CdWpTtChi2RecTool::finalize() {
    if (!m_reg.finalize()) return false;
    LogInfo << m_name << " finalized successfully\n";
    return true;
}