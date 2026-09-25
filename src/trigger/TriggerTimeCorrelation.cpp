#include "trigger/TriggerTimeCorrelation.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolBase.h"

void TriggerTimeCorrelation::configure(const SniperJSON& config) {
    if (!config.valid()) return;
    setConfigValue(m_time_window.first, "LowerBound", config);
    setConfigValue(m_time_window.second, "UpperBound", config);
}

bool TriggerTimeCorrelation::correct(RecPmtTable& table, JM::NavBuffer* buf) {
    LogDebug << "Navigation buffer: " << buf << '\n';
    (void)table;
    std::cin.get();
    return true;
}