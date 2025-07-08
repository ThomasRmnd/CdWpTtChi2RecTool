#include "transformer/PmtLocTypeTransformer.h"

#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(PmtLocTypeTransformer);

void PmtLocTypeTransformer::transform(RecPmtProp& pmt) {
    if (!pmt.used) return;
    pmt.used = !checkPmtType(pmt);
    return;
}