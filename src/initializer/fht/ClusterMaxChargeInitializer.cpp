#include "initializer/fht/ClusterMaxChargeInitializer.hpp"

#include <unordered_map>

#include <TH1D.h>

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(ClusterMaxChargeInitializer);

ClusterMaxChargeInitializer::ClusterMaxChargeInitializer(const std::string& name) : 
    Initializer<FhtMethodTag>{name},
    m_pmt_cnt_thold{5.0},
    m_shift{2.0},
    m_range{5.0},
    m_pe_thold{10.0},
    m_radius{9418.0},
    m_lwr_q_thold{0.9},
    m_dir_corr_factor{1.5},
    m_radius_end{10000.0},
    m_hist{std::make_unique<TH1D>("h__ClusterMaxChargeInitializer", "h__ClusterMaxChargeInitializer", 1000, 0.0, 1000.0)}
{
    m_hist->SetDirectory(0);
}

ClusterMaxChargeInitializer::ClusterMaxChargeInitializer(const std::string& name, double pmt_cnt_thold, double shift, double range, double pe_thold, double radius, double lwr_q_thold, double dir_corr_factor, double radius_end) : 
    Initializer<FhtMethodTag>{name},
    m_pmt_cnt_thold{pmt_cnt_thold},
    m_shift{shift},
    m_range{range},
    m_pe_thold{pe_thold},
    m_radius{radius},
    m_lwr_q_thold{lwr_q_thold},
    m_dir_corr_factor{dir_corr_factor},
    m_radius_end{radius_end},
    m_hist{std::make_unique<TH1D>("h__ClusterMaxChargeInitializer", "h__ClusterMaxChargeInitializer", 1000, 0.0, 1000.0)}
{
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
    while (m_hist->GetBinContent(idx) > m_pmt_cnt_thold && idx > 0) {
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
        if (m_radius_end < mag(it->pos - pos_q_cntrd) || mag(it->pos - ipos) < m_radius) continue;
        if (max_q < it->q) max_q = it->q;
    }
    vec3 fpos;
    std::size_t n = 0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (m_radius_end < mag(it->pos - pos_q_cntrd) || mag(it->pos - ipos) < m_radius) continue;
        if (it->q < m_lwr_q_thold * max_q) continue;
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