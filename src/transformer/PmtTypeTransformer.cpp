#include "transformer/PmtTypeTransformer.hpp"

PmtTypeTransformer::PmtTypeTransformer(const std::string& name, const RecPmtType& pmt_type) :
    Transformer{name, pmt_type}
{}

void PmtTypeTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used) return;
    pmt.used = !checkPmtType(pmt);
}