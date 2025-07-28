#include "transformer/PmtLocTypeTransformer.hpp"

PmtLocTypeTransformer::PmtLocTypeTransformer(const std::string& name, const RecPmtType& pmt_type) :
    Transformer{name, pmt_type}
{}

void PmtLocTypeTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used) return;
    pmt.used = !checkPmtType(pmt);
    return;
}