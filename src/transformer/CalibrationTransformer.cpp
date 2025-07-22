#include "transformer/CalibrationTransformer.hpp"

#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(CalibrationTransformer);

CalibrationTransformer::CalibrationTransformer(const std::string& name) :
    Transformer(name),
    m_offset(0.0)
{}

CalibrationTransformer::CalibrationTransformer(const std::string& name, const RecPmtType& type, double offset) :
    Transformer(name, type),
    m_offset(offset)
{}

void CalibrationTransformer::configure(const SniperJSON& config) {
    LogDebug << "Before: TimeOffset = " << m_offset << '\n';
    m_offset = getConfigValue<double>("TimeOffset", config);
    LogDebug << "After: TimeOffset = " << m_offset << '\n';
}

void CalibrationTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    pmt.fht += m_offset;
    return;
}