#include "transformer/PmtLocTypeTransformer.hpp"

#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(PmtLocTypeTransformer);

void PmtLocTypeTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used) return;
    pmt.used = !checkPmtType(pmt);
    return;
}