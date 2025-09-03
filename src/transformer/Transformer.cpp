#include "transformer/Transformer.hpp"

#include "SniperKernel/SniperLog.h"

Transformer::Transformer(const std::string& name, const RecPmtType& pmt_type) :
    Configurable{name},
    PmtTypeChecker{pmt_type}
{}

void Transformer::transform(RecPmtTable& table) {
    if (!findRange(table)) return;

    std::size_t isize = table.size();
    table.erase(std::remove_if(m_ftable, m_ltable, [&](RecPmtProp& pmt) { transformPmt(pmt); return !pmt.used; }), m_ltable);
    std::size_t fsize = table.size();

    LogDebug << isize << " -> " << fsize << " = " << isize - fsize << " PMTs are removed\n";
}

bool Transformer::findRange(RecPmtTable& table) {
    m_ftable = std::find_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return checkPmtType(pmt); });
    m_ltable = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return checkPmtType(pmt); }).base();

    if (std::distance(m_ftable, m_ltable) <= 0) {
        LogDebug << "No PMTs match the type for this Transformer\n";
        return false;
    }
    return true;
}