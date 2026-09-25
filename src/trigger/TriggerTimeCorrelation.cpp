#include "trigger/TriggerTimeCorrelation.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolBase.h"

TriggerTimeCorrelation::TriggerTimeCorrelation(const std::string& name, const std::pair<double, double>& time_window) :
    Configurable{name},
    m_time_window{time_window}
{}

void TriggerTimeCorrelation::configure(const SniperJSON& config) {
    if (!config.valid()) return;
    setConfigValue(m_time_window.first, "LowerBound", config);
    setConfigValue(m_time_window.second, "UpperBound", config);
}

bool TriggerTimeCorrelation::correct(RecPmtTable& table, JM::NavBuffer* buf) {
    if (!buf) {
        LogError << "NavBuffer is nullptr\n";
        return false;
    }
    JM::EvtNavigator* nav = buf->curEvt();
    if (!nav) {
        LogError << "Current EvtNavigator is nullptr\n";
        return false;
    }

    m_cur_ts = TimeStamp{nav->TimeStamp().GetTimeSpec()};
    m_cur_type = getDetectorType(nav);
    if (m_cur_type == DetectorType::NONE) {
        LogInfo << "Unexpected entry type (" << nav->getDetectorType() << ") at " << m_cur_ts << ". Ignoring\n";
        return false;
    }

    for (JM::NavBuffer::Iterator it = buf->begin(); it != buf->end(); ++it) {
        if (it == buf->current()) continue;
        nav = it->get();
        if (!nav) {
            LogError << "EvtNavigator is nullptr\n";
            return false;
        }
        m_other_ts = TimeStamp{nav->TimeStamp().GetTimeSpec()};
        m_other_type = getDetectorType(nav);
        if (m_other_type == DetectorType::NONE) {
            LogInfo << "Unexpected entry type (" << nav->getDetectorType() << ") at " << m_other_ts << ". Ignoring\n";
            continue;
        }
        LogDebug << "Current entry: (type: " << static_cast<int>(m_cur_type) << ", ts: " << m_cur_ts << ")\n";
        LogDebug << "Other entry: (type: " << static_cast<int>(m_other_type) << ", ts: " << m_other_ts << ")\n";
        if ( (m_cur_type & m_other_type) != DetectorType::NONE ) continue; // current event has already assotiated an event with same type
        TimeStamp diff_ts = m_cur_ts - m_other_ts;
        double diff_ts_ns = static_cast<double>(diff_ts.GetSec()) * 1.0e9 + static_cast<double>(diff_ts.GetNanoSec());
        if (diff_ts_ns < m_time_window.first || m_time_window.second < diff_ts_ns) continue;
        LogInfo << "Entries are associated: " << m_cur_ts << " - " << m_other_ts << " = " << diff_ts_ns << " ns\n";

        changeRefTime(table, diff_ts_ns);
        m_cur_type |= m_other_type;
    }

    return true;
}

DetectorType TriggerTimeCorrelation::getDetectorType(JM::EvtNavigator* nav) {
    DetectorType type = DetectorType::NONE;

    JM::EvtNavigator::DetectorType evt_type = nav->getDetectorType();

    if (evt_type == JM::EvtNavigator::DetectorType::CD) type |= DetectorType::CD;
    if (evt_type == JM::EvtNavigator::DetectorType::WP) type |= DetectorType::WP;
    if (evt_type == JM::EvtNavigator::DetectorType::TT) type |= DetectorType::TT;

    return type;
}

void TriggerTimeCorrelation::changeRefTime(RecPmtTable& table, double diff_ts) {
    if ( (m_other_type & DetectorType::CD) == DetectorType::CD ) changeRefTimeInRange(table, diff_ts, RecPmtType::PMT_CD);
    if ( (m_other_type & DetectorType::WP) == DetectorType::WP ) changeRefTimeInRange(table, diff_ts, RecPmtType::PMT_WP);
    if ( (m_other_type & DetectorType::TT) == DetectorType::TT ) changeRefTimeInRange(table, diff_ts, RecPmtType::PMT_TT);
}

void TriggerTimeCorrelation::changeRefTimeInRange(RecPmtTable& table, double diff_ts, const RecPmtType& type) {
    for (RecPmtProp& pmt : table) {
        if (!pmt.used || !hasPmtType(pmt, type)) continue;
        LogDebug << "PMT (" << pmt.id << " = " << pmt.fht << " --> " << pmt.fht - diff_ts << '\n';
        pmt.fht -= diff_ts;
    }
}