#include "transformer/CalibTimeTransformer.hpp"

CalibTimeTransformer::CalibTimeTransformer(const std::string& name, const RecPmtType& type, double offset) :
    Transformer(name, type),
    m_offset(offset)
{}

void CalibTimeTransformer::configure(const SniperJSON& config) {
    setConfigValue(m_offset, "TimeOffset", config);
}

void CalibTimeTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    pmt.fht += m_offset;
    return;
}