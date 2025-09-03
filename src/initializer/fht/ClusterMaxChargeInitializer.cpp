#include "initializer/fht/ClusterMaxChargeInitializer.hpp"

#include <unordered_map>

#include <TH1D.h>

#include "SniperKernel/SniperLog.h"

ClusterMaxChargeInitializer::ClusterMaxChargeInitializer(const std::string& name, double pmt_cnt_thold, double shift, double range, double pe_thold, double ipos_radius, double q_ratio, double dir_corr_factor, double fpos_radius) : 
    Initializer<FhtMethodTag>{name},
    m_pmt_cnt_thold{pmt_cnt_thold},
    m_shift{shift},
    m_range{range},
    m_pe_thold{pe_thold},
    m_ipos_radius{ipos_radius},
    m_q_ratio{q_ratio},
    m_dir_corr_factor{dir_corr_factor},
    m_fpos_radius{fpos_radius},
    m_hist{std::make_unique<TH1D>((m_name + "__Histogram").c_str(), (m_name + "__Histogram").c_str(), 1000, 0.0, 1000.0)}
{
    m_hist->SetDirectory(0);
}

void ClusterMaxChargeInitializer::configure(const SniperJSON& config) {
    if (!config.valid()) return;
    setConfigValue(m_pmt_cnt_thold, "PmtThreshold", config);
    setConfigValue(m_shift, "TimeShift", config);
    setConfigValue(m_range, "TimeRange", config);
    setConfigValue(m_pe_thold, "ChargeThreshold", config);
    setConfigValue(m_ipos_radius, "IposRadius", config);
    setConfigValue(m_q_ratio, "ChargeRatio", config);
    setConfigValue(m_dir_corr_factor, "DirectionCorrectionFactor", config);
    setConfigValue(m_fpos_radius, "FposRadius", config);
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

bool ClusterMaxChargeInitializer::getTable(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) {
    m_it_table.clear();
    m_hist->Reset();

    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        m_hist->Fill(it->fht);
    }

    int idx = m_hist->GetMaximumBin();
    while (m_hist->GetBinContent(idx) > m_pmt_cnt_thold && idx > 1) {
        --idx;
    }
    double itime = m_hist->GetBinCenter(idx) - m_shift;

    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (m_range < std::abs(it->fht - itime)) continue;
        if (it->q < m_pe_thold) continue;
        m_it_table.push_back(it);
    }
    if (m_it_table.empty()) {
        LogWarn << "No valid PMTs found in the time window\n";
        return false;
    }
    return true;
}

double ClusterMaxChargeInitializer::getITime() {
    double itime = 0.0;
    std::size_t n = 0ul;
    for (RecPmtTable::const_iterator it : m_it_table) {
        itime += it->fht;
        ++n;
    }
    return itime / static_cast<double>(n);
}

vec3 ClusterMaxChargeInitializer::getIPos() {
    vec3 ipos;
    std::size_t n = 0ul;
    for (RecPmtTable::const_iterator it : m_it_table) {
        ipos += it->pos;
        ++n;
    }
    return ipos / static_cast<double>(n);
}

vec3 ClusterMaxChargeInitializer::getChargeCentroid(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) {
    vec3 pos;
    double totq = 0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        pos += it->pos * it->totq;
        totq += it->totq;
    }
    return pos / totq;
}

vec3 ClusterMaxChargeInitializer::getFPos(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable, const vec3& ipos) {
    vec3 pos_q_cntrd = m_dir_corr_factor * getChargeCentroid(ftable, ltable);
    vec3 dir = unit(pos_q_cntrd - ipos);
    pos_q_cntrd = ipos - 2.0 * dot(ipos, dir) * dir;
    double max_q = 0.0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (m_fpos_radius < mag(it->pos - pos_q_cntrd) || mag(it->pos - ipos) < m_ipos_radius) continue;
        if (max_q < it->q) max_q = it->q;
    }
    vec3 fpos;
    std::size_t n = 0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (m_fpos_radius < mag(it->pos - pos_q_cntrd) || mag(it->pos - ipos) < m_ipos_radius) continue;
        if (it->q < m_q_ratio * max_q) continue;
        fpos += it->pos;
        ++n;
    }
    if (!n) {
        LogWarn << "No valid PMTs found for the end point\n";
        return pos_q_cntrd;
    }
    return fpos / static_cast<double>(n);
}

bool ClusterMaxChargeInitializer::initiate(const RecPmtTable& table) {
    unsigned int count = std::count_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_20INCH); });
    LogDebug << count << " PMTs are used for the initialization\n";

    RecPmtTable::const_iterator ftable = std::find_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_20INCH); });
    RecPmtTable::const_iterator ltable = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_20INCH); }).base();

    if (!getTable(ftable, ltable)) return false;
    double itime = getITime();
    vec3 ipos = getIPos();
    vec3 fpos = getFPos(ftable, ltable, ipos);
    vec3 dir = fpos - ipos;
    // Preliminary
    vec3 mid = ipos - dot(ipos, dir) * dir / mag2(dir);
    mid /= 1.0225; // = (1.025 + 1.020) / 2.0 // best guess for now --> might need to improve init point for better angle 
    dir = mid - ipos;
    m_params = std::vector<double>{itime, theta(ipos), phi(ipos), theta(dir), phi(dir)};

    return true;
}

ParamsType ClusterMaxChargeInitializer::getOParamsType() {
    return ParamsType::SingleAcrylic;
}