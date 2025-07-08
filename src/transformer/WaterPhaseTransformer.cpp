#include "transformer/WaterPhaseTransformer.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(WaterPhaseTransformer);

#ifdef __USE_WATERPHASETRANSFORMER_J24_3_0_SIMULATIONS_IMPLEMENTATION__

WaterPhaseTransformer::WaterPhaseTransformer(const std::string& name, double q_thold, int nb_bins, double xmin, double xmax, double lpmt_thold_min, double lpmt_thold_max, double spmt_thold_min, double spmt_thold_max, double radius, unsigned int nb_neigh_thold, double neigh_fht_diff) :
    Transformer(name, RecPmtType::PMT_CD),
    m_hist(std::make_unique<TH1D>("fht_hist", "FHT Distribution", nb_bins, xmin, xmax)),
    m_q_thold(q_thold),
    m_lpmt_thold_min(lpmt_thold_min),
    m_lpmt_thold_max(lpmt_thold_max),
    m_spmt_thold_min(spmt_thold_min),
    m_spmt_thold_max(spmt_thold_max),
    m_r2(radius * radius),
    m_nb_neigh_thold(nb_neigh_thold),
    m_neigh_fht_diff(neigh_fht_diff)
{
    m_hist->SetDirectory(0);
}

void WaterPhaseTransformer::getTimes20inch(const RecPmtTable& table) {
    m_hist->Reset();
    for (const RecPmtProp& pmt : table) {
        if (!pmt.used || (pmt.type & RecPmtType::PMT_20INCH) != pmt.type) continue;
        m_hist->Fill(pmt.fht, std::log10(pmt.hitq[0]));
    }

    int idx = m_hist->GetMaximumBin();
    while (m_lpmt_thold_min < m_hist->GetBinContent(idx) && 1 < idx) {
        --idx;
    }
    m_lpmt_t_i = m_hist->GetBinCenter(idx);

    idx = m_hist->GetMaximumBin();
    while (m_lpmt_thold_max < m_hist->GetBinContent(idx) && idx < m_hist->GetNbinsX()) {
        ++idx;
    }
    m_lpmt_t_f = m_hist->GetBinCenter(idx);
}

void WaterPhaseTransformer::getTimes3inch(const RecPmtTable& table) {
    m_hist->Reset();
    for (const RecPmtProp& pmt : table) {
        if (!pmt.used || (pmt.type & RecPmtType::PMT_3INCH) != pmt.type) continue;
        m_hist->Fill(pmt.fht);
    }

    int idx = m_hist->GetMaximumBin();
    while (m_spmt_thold_min < m_hist->GetBinContent(idx) && 1 < idx) {
        --idx;
    }
    m_spmt_t_i = m_hist->GetBinCenter(idx);

    idx = m_hist->GetMaximumBin();
    while (m_spmt_thold_max < m_hist->GetBinContent(idx) && idx < m_hist->GetNbinsX()) {
        ++idx;
    }
    m_spmt_t_f = m_hist->GetBinCenter(idx);
}

void WaterPhaseTransformer::operator()(RecPmtTable& table) {
    if (!findRange(table)) return;
    
    m_nb_spmt = 0u;
    m_nb_lpmt = 0u;

    std::size_t isize_lpmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_20INCH) == pmt.type; });
    std::size_t isize_spmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });

    for (RecPmtProp& pmt : table) {
        if (!pmt.used || (pmt.type & RecPmtType::PMT_20INCH) != pmt.type) continue;
        if (pmt.hitq[0] < m_q_thold || 10000.0 < pmt.q) pmt.used = false; // 10000.0 because of slight bug in the simulation
    }
    getTimes20inch(table);
    getTimes3inch(table);

    std::for_each(m_ftable, m_ltable, [&](RecPmtProp& pmt) { transform(pmt); });

    std::size_t tsize_spmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });
    std::cout << "SPMT: " << isize_spmt << " - " << tsize_spmt << " = " << isize_spmt - tsize_spmt << '\n';
    std::cout << m_nb_spmt << " SPMTs and " << m_nb_lpmt << " LPMTs are used\n";

    // Neighbouring filter for 3-inch PMTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m_nb_neigh.clear();
    m_nb_neigh.resize(m_nb_spmt, 0u);
    m_fht_neigh.clear();
    m_fht_neigh.resize(m_nb_spmt, 0.0);

    std::size_t k = 0u;
    for (RecPmtProp& ipmt : table) {
        if (!ipmt.used || (ipmt.type & RecPmtType::PMT_3INCH) != ipmt.type) continue;
        for (RecPmtProp& jpmt : table) {
            if (!jpmt.used || (jpmt.type & RecPmtType::PMT_3INCH) != jpmt.type) continue;
            if (ipmt.id == jpmt.id) continue;
            if (m_r2 < mag2(ipmt.pos - jpmt.pos)) continue;
            m_nb_neigh[k] += 1u;
            m_fht_neigh[k] += jpmt.fht;
        }
        if (m_nb_neigh[k] != 0u) m_fht_neigh[k] /= static_cast<double>(m_nb_neigh[k]);
        ++k;
    }

    std::cout << k << '\n';

    k = 0u;
    for (RecPmtProp& pmt : table) {
        if (!pmt.used || (pmt.type & RecPmtType::PMT_3INCH) != pmt.type) continue;
        if (m_nb_neigh[k] < m_nb_neigh_thold) pmt.used = false;
        if ( m_neigh_fht_diff < std::abs(pmt.fht - m_fht_neigh[k]) ) pmt.used = false;
        ++k;
    }

    std::cout << k << '\n';

    // Neighbouring filter for 20-inch PMTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m_nb_neigh.clear();
    m_nb_neigh.resize(m_nb_lpmt, 0u);
    m_fht_neigh.clear();
    m_fht_neigh.resize(m_nb_lpmt, 0.0);
    double min_fht = std::numeric_limits<double>::infinity();
    k = 0u;
    for (RecPmtProp& ipmt : table) {
        if (!ipmt.used || (ipmt.type & RecPmtType::PMT_20INCH) != ipmt.type) continue;
        for (RecPmtProp& jpmt : table) {
            if (!jpmt.used || (jpmt.type & RecPmtType::PMT_20INCH) != jpmt.type) continue;
            if (ipmt.id == jpmt.id) continue;
            if (m_r2 < mag2(ipmt.pos - jpmt.pos)) continue;
            m_nb_neigh[k] += 1u;
            m_fht_neigh[k] += jpmt.fht;
        }
        if (ipmt.fht < min_fht) min_fht = ipmt.fht;
        if (m_nb_neigh[k] != 0u) m_fht_neigh[k] /= static_cast<double>(m_nb_neigh[k]);
        ++k;
    }

    k = 0u;
    for (RecPmtProp& pmt : table) {
        if (!pmt.used || (pmt.type & RecPmtType::PMT_20INCH) != pmt.type) continue;
        if (m_nb_neigh[k] < 1 && min_fht + 50.0 < pmt.fht) pmt.used = false;
        ++k;
    }

    std::size_t isize = table.size();
    table.erase(std::remove_if(m_ftable, m_ltable, [&](RecPmtProp& pmt) { return !pmt.used; }), m_ltable);
    std::size_t fsize = table.size();

    std::size_t fsize_lpmt = std::count_if(table.begin(), table.end(), [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_20INCH) == pmt.type; });
    std::cout << "LPMT: " << isize_lpmt << " - " << fsize_lpmt << " = " << isize_lpmt - fsize_lpmt << '\n';

    std::size_t fsize_spmt = std::count_if(table.begin(), table.end(), [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });
    std::cout << "SPMT: " << tsize_spmt << " - " << fsize_spmt << " = " << tsize_spmt - fsize_spmt << '\n';
        
    LogDebug << isize << " -> " << fsize << " = " << isize - fsize << " PMTs are removed." << '\n';

    return;
}

void WaterPhaseTransformer::transform(RecPmtProp& pmt) {
    if (!pmt.used || !check(pmt)) return;
    if ( (pmt.type & RecPmtType::PMT_20INCH) == pmt.type ) {
        if ( (pmt.fht < m_lpmt_t_i && pmt.hitq[0] < m_q_thold) || m_lpmt_t_f < pmt.fht ) pmt.used = false;
        else ++m_nb_lpmt;
    }
    else if ( (pmt.type & RecPmtType::PMT_3INCH) == pmt.type ) {
        if ( pmt.fht < m_spmt_t_i || m_spmt_t_f < pmt.fht ) pmt.used = false;
        else ++m_nb_spmt;
    }
    return;
}

#else

WaterPhaseTransformer::WaterPhaseTransformer(const std::string& name, int nb_bins, double xmin, double xmax, double q_thold_itime, double ratio_entries_ftime, double r, unsigned int neigh_thold, double dt) :
    Transformer(name, RecPmtType::PMT_CD),
    m_hist(std::make_unique<TH1D>("fht_hist", "FHT Distribution", nb_bins, xmin, xmax)),
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
        if (!pmt.used) continue;
        if ( (pmt.type & RecPmtType::PMT_20INCH) != pmt.type ) continue;

#if __USE_RECPMTPROP_VERSION__ == 1
        if (static_cast<double>(pmt.hitq[0]) < m_q_thold_itime || pmt.fht < 100.0) continue;
#elif __USE_RECPMTPROP_VERSION__ == 2
        if (static_cast<double>(pmt.it_tq->q) < m_q_thold_itime || pmt.fht < 100.0) continue;
#elif __USE_RECPMTPROP_VERSION__ == 3
        if (pmt.q < m_q_thold_itime || pmt.fht < 100.0) continue;
#endif // __USE_RECPMTPROP_VERSION__

        m_lpmt_itime = std::min(m_lpmt_itime, pmt.fht);
        flag = true;
    }
    if (!flag) LogWarn << "No valid PMTs found to get the initial time\n";
    return flag;
}

bool WaterPhaseTransformer::getFTime20inch(const RecPmtTable& table) {
    m_hist->Reset();
    for (const RecPmtProp& pmt : table) {
        if (!pmt.used) continue;
        if ( (pmt.type & RecPmtType::PMT_20INCH) != pmt.type ) continue;
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

void WaterPhaseTransformer::operator()(RecPmtTable& table) {
    if (!findRange(table)) return;

    std::size_t isize_lpmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_20INCH) == pmt.type; });
    std::size_t isize_spmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });

    if (!getITime20inch(table)) return;
    if (!getFTime20inch(table)) return;

    std::cout << "itime = " << m_lpmt_itime << ", ftime = " << m_lpmt_ftime << '\n';

    std::for_each(m_ftable, m_ltable, [&](RecPmtProp& pmt) { transform(pmt); });

    std::size_t tsize_lpmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_20INCH) == pmt.type; });
    std::size_t tsize_spmt = std::count_if(m_ftable, m_ltable, [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });
    
    std::cout << "LPMT: " << isize_lpmt << " - " << tsize_lpmt << " = " << isize_lpmt - tsize_lpmt << '\n';
    std::cout << "SPMT: " << isize_spmt << " - " << tsize_spmt << " = " << isize_spmt - tsize_spmt << '\n';

    // Neighbouring filter for 20-inch PMTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m_nb_neigh.clear();
    m_nb_neigh.resize(table.size(), 0u);

    for (std::size_t i = 0; i < table.size(); ++i) {
        if (!table[i].used) continue;
        if ( (table[i].type & RecPmtType::PMT_20INCH) != table[i].type ) continue;
        for (std::size_t j = 0; j < table.size(); ++j) {
            if (!table[j].used) continue;
            if ( (table[j].type & RecPmtType::PMT_20INCH) != table[j].type ) continue;
            if (i == j) continue;
            if (m_r2 < mag2(table[i].pos - table[j].pos)) continue;
            ++m_nb_neigh[i];
        }
    }

    for (std::size_t k = 0; k < table.size(); ++k) {
        if (!table[k].used) continue;
        if ( (table[k].type & RecPmtType::PMT_20INCH) != table[k].type ) continue;
        if (m_nb_neigh[k] < m_neigh_thold && m_lpmt_itime + m_dt < table[k].fht) table[k].used = false;
    }

    std::size_t isize = table.size();
    table.erase(std::remove_if(m_ftable, m_ltable, [&](RecPmtProp& pmt) { return !pmt.used; }), m_ltable);
    std::size_t fsize = table.size();

    std::size_t fsize_lpmt = std::count_if(table.begin(), table.end(), [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_20INCH) == pmt.type; });
    std::size_t fsize_spmt = std::count_if(table.begin(), table.end(), [](const RecPmtProp& pmt) { return pmt.used && (pmt.type & RecPmtType::PMT_3INCH) == pmt.type; });
    
    std::cout << "LPMT: " << isize_lpmt << " - " << fsize_lpmt << " = " << isize_lpmt - fsize_lpmt << '\n';
    std::cout << "SPMT: " << tsize_spmt << " - " << fsize_spmt << " = " << tsize_spmt - fsize_spmt << '\n';
        
    LogDebug << isize << " -> " << fsize << " = " << isize - fsize << " PMTs are removed\n";

    return;
}

void WaterPhaseTransformer::transform(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    if ( (pmt.type & RecPmtType::PMT_20INCH) == pmt.type ) {
        if ( pmt.fht < m_lpmt_itime || m_lpmt_ftime < pmt.fht ) pmt.used = false;
    }
    // else if ( (pmt.type & RecPmtType::PMT_3INCH) == pmt.type ) {}
    return;
}

#endif // __USE_WATERPHASETRANSFORMER_J24_3_0_SIMULATIONS_IMPLEMENTATION__