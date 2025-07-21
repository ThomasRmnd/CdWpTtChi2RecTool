#include "transformer/FhtChargeTholdTransformer.hpp"

#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(FhtChargeTholdTransformer);

FhtChargeTholdTransformer::FhtChargeTholdTransformer(const std::string& name) :
    Transformer(name),
    m_q_thold(0.0)
{}

FhtChargeTholdTransformer::FhtChargeTholdTransformer(const std::string& name, const RecPmtType& type, double q_thold) :
    Transformer(name, type),
    m_q_thold(q_thold)
{}

void FhtChargeTholdTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    pmt.used = (pmt.q >= m_q_thold);
    return;
}