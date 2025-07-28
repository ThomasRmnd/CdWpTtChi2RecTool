#include "transformer/CalibrationTransformer.hpp"

DECLARE_TOOL(CalibrationTransformer);

CalibrationTransformer::CalibrationTransformer(const std::string& name) :
    Transformer(name, RecPmtType::PMT_UNKNOWN),
    m_offset(0.0)
{}

CalibrationTransformer::CalibrationTransformer(const std::string& name, const RecPmtType& type, double offset) :
    Transformer(name, type),
    m_offset(offset)
{}

void CalibrationTransformer::configure(const SniperJSON& config) {
    setConfigValue(m_offset, "TimeOffset", config);
}

void CalibrationTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    pmt.fht += m_offset;
    return;
}