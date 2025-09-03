#include "transformer/EarlyLateFhtTransformer.hpp"

#include "SniperKernel/SniperLog.h"

EarlyLateFhtTransformer::EarlyLateFhtTransformer(const std::string& name, const RecPmtType& type, int nb_bins, double xmin, double xmax, double pmt_thold, double shift, double relative_cut) :
    Transformer(name, type),
    m_hist(std::make_unique<TH1D>((m_name + "__Histogram").c_str(), (m_name + "__Histogram").c_str(), nb_bins, xmin, xmax)),
    m_pmt_thold(pmt_thold),
    m_shift(shift),
    m_relative_cut(relative_cut),
    m_itime(0.0)
{
    m_hist->SetDirectory(0);
}

void EarlyLateFhtTransformer::configure(const SniperJSON& config) {
    if (!config.valid()) return;
    setConfigValue(m_pmt_thold, "PmtThreshold", config);
    setConfigValue(m_shift, "TimeShift", config);
    setConfigValue(m_relative_cut, "RelativeTimeCut", config);
    int nbins = 0;
    double xmin = 0.0, xmax = 0.0;
    if (
        !setConfigValue(nbins, "Histogram_Nbins", config) ||
        !setConfigValue(xmin, "Histogram_Xmin", config) ||
        !setConfigValue(xmax, "Histogram_Xmax", config)
    ) return;
    m_hist = std::make_unique<TH1D>((m_name + "__Histogram").c_str(), (m_name + "__Histogram").c_str(), nbins, xmin, xmax);
    m_hist->SetDirectory(0);
}

double EarlyLateFhtTransformer::getITime(const RecPmtTable& table) {
    m_hist->Reset();
    for (const RecPmtProp& pmt : table) {
        if (!pmt.used || !checkPmtType(pmt)) continue;
        m_hist->Fill(pmt.fht);
    }
    
    int idx = m_hist->GetMaximumBin();
    while (m_pmt_thold < m_hist->GetBinContent(idx) && 0 < idx) {
        --idx;
    }

    return (m_hist->GetBinCenter(idx) - m_shift);
}

void EarlyLateFhtTransformer::transform(RecPmtTable& table) {
    m_itime = getITime(table);
    Transformer::transform(table);
    return;
}

void EarlyLateFhtTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    pmt.used = (m_itime <= pmt.fht && pmt.fht <= m_itime + m_relative_cut);
    return;
}