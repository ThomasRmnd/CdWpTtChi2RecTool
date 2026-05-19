#include "transformer/WaterPhaseTransformer.hpp"

#include "SniperKernel/SniperLog.h"

WaterPhaseTransformer::WaterPhaseTransformer(const std::string& name, int nb_bins, double xmin, double xmax, double q_thold_itime, double ratio_entries_ftime, double r, unsigned int neigh_thold, double dt) :
    Transformer(name, RecPmtType::PMT_20INCH),
    m_hist(std::make_unique<TH1D>("h__WaterPhaseTransformer", "h__WaterPhaseTransformer", nb_bins, xmin, xmax)),
    m_q_thold_itime(q_thold_itime),
    m_ratio_entries_ftime(ratio_entries_ftime),
    m_r2(r * r),
    m_neigh_thold(neigh_thold),
    m_dt(dt)
{
    m_hist->SetDirectory(0);
}

bool WaterPhaseTransformer::getITime20inch(const RecPmtTable& table) {
    m_lpmt_itime = std::numeric_limits<double>::infinity();
    bool flag = false;
    for (const RecPmtProp& pmt : table) {
        if (!pmt.used || !checkPmtType(pmt)) continue;
        if (pmt.q < m_q_thold_itime || pmt.fht < 100.0) continue;
        m_lpmt_itime = std::min(m_lpmt_itime, pmt.fht);
        flag = true;
    }
    if (!flag) LogWarn << "No valid PMTs found to get the initial time\n";
    return flag;
}

bool WaterPhaseTransformer::getFTime20inch(const RecPmtTable& table) {
    m_hist->Reset();
    for (const RecPmtProp& pmt : table) {
        if (!pmt.used || !checkPmtType(pmt)) continue;
        m_hist->Fill(pmt.fht);
    }

    if (m_hist->GetEntries() == 0) {
        LogWarn << "No valid PMTs found to get the final time\n";
        return false;
    }

    int idx = m_hist->GetMaximumBin();
    double entries_thold = m_ratio_entries_ftime * m_hist->GetMaximum();
    while (entries_thold < m_hist->GetBinContent(idx) && idx < m_hist->GetNbinsX()) {
        ++idx;
    }
    m_lpmt_ftime = m_hist->GetBinCenter(idx);
    return true;
}

void WaterPhaseTransformer::transform(RecPmtTable& table) {
    if (!findRange(table)) return;

    std::size_t isize_lpmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_20INCH) == pmt.type; });
    std::size_t isize_spmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });

    if (!getITime20inch(table)) return;
    if (!getFTime20inch(table)) return;

    LogDebug << "itime = " << m_lpmt_itime << ", ftime = " << m_lpmt_ftime << '\n';

    std::for_each(m_ftable, m_ltable, [&](RecPmtProp& pmt) { transformPmt(pmt); });

    std::size_t tsize_lpmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_20INCH) == pmt.type; });
    std::size_t tsize_spmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });
    
    LogDebug << "LPMT: " << isize_lpmt << " - " << tsize_lpmt << " = " << isize_lpmt - tsize_lpmt << '\n';
    LogDebug << "SPMT: " << isize_spmt << " - " << tsize_spmt << " = " << isize_spmt - tsize_spmt << '\n';

    // Neighbouring filter for 20-inch PMTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m_nb_neigh.clear();
    m_nb_neigh.resize(table.size(), 0u);

    for (std::size_t i = 0; i < table.size(); ++i) {
        if (!table[i].used || !checkPmtType(table[i])) continue;
        for (std::size_t j = 0; j < table.size(); ++j) {
            if (!table[j].used || !checkPmtType(table[j])) continue;
            if (i == j) continue;
            if (m_r2 < mag2(table[i].pos - table[j].pos)) continue;
            ++m_nb_neigh[i];
        }
    }

    for (std::size_t k = 0; k < table.size(); ++k) {
        if (!table[k].used || !checkPmtType(table[k])) continue;
        if (m_nb_neigh[k] < m_neigh_thold && m_lpmt_itime + m_dt < table[k].fht) table[k].used = false;
    }

    std::size_t isize = table.size();
    table.erase(std::remove_if(m_ftable, m_ltable, [&](RecPmtProp& pmt) { return !pmt.used; }), m_ltable);
    std::size_t fsize = table.size();

    std::size_t fsize_lpmt = std::count_if(table.begin(), table.end(), [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_20INCH) == pmt.type; });
    std::size_t fsize_spmt = std::count_if(table.begin(), table.end(), [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });
    
    LogDebug << "LPMT: " << isize_lpmt << " - " << fsize_lpmt << " = " << isize_lpmt - fsize_lpmt << '\n';
    LogDebug << "SPMT: " << tsize_spmt << " - " << fsize_spmt << " = " << tsize_spmt - fsize_spmt << '\n';
        
    LogDebug << isize << " -> " << fsize << " = " << isize - fsize << " PMTs are removed\n";

    return;
}

void WaterPhaseTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    if ( pmt.fht < m_lpmt_itime || m_lpmt_ftime < pmt.fht ) pmt.used = false;
    return;
}