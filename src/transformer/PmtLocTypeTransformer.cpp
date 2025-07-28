#include "transformer/PmtLocTypeTransformer.hpp"

#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(PmtLocTypeTransformer);

PmtLocTypeTransformer::PmtLocTypeTransformer(const std::string& name) :
    Transformer{name, RecPmtType::PMT_UNKNOWN}
{}

PmtLocTypeTransformer::PmtLocTypeTransformer(const std::string& name, const RecPmtType& pmt_type) :
    Transformer{name, pmt_type}
{}

void PmtLocTypeTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used) return;
    pmt.used = !checkPmtType(pmt);
    return;
}