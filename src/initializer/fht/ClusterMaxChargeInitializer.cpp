#include "initializer/fht/ClusterMaxChargeInitializer.hpp"

#include <unordered_map>

#include <TH1D.h>

#include "SniperKernel/SniperLog.h"

ClusterMaxChargeInitializer::ClusterMaxChargeInitializer(const std::string& name, int nbins_cd, int nbins_wp, double xmin, double xmax, double window_min, double window_max, double qthold, double offset_final_cd, double exclusion_time_wp, double leading_entries_ratio, Mode mode) :
    Initializer<FhtMethodTag>(name),
    m_hist_cd(std::make_unique<TH1D>((m_name + "__Histogram_CD").c_str(), (m_name + "__Histogram_CD").c_str(), nbins_cd, xmin, xmax)),
    m_hist_wp(std::make_unique<TH1D>((m_name + "__Histogram_WP").c_str(), (m_name + "__Histogram_WP").c_str(), nbins_wp, xmin, xmax)),
    m_window_min(window_min),
    m_window_max(window_max),
    m_qthold(qthold),
    m_offset_final_cd(offset_final_cd),
    m_exclusion_time_wp(exclusion_time_wp),
    m_leading_entries_ratio(leading_entries_ratio),
    m_mode(mode)
{
    m_hist_cd->SetDirectory(0);
    m_hist_wp->SetDirectory(0);
}

void ClusterMaxChargeInitializer::configure(const SniperJSON& config) {
    if (!config.valid()) return;
    setConfigValue(m_window_min, "WindowMin", config);
    setConfigValue(m_window_max, "WindowMax", config);
    setConfigValue(m_qthold, "ChargeThreshold", config);
    setConfigValue(m_offset_final_cd, "CDFinalOffset", config);
    setConfigValue(m_exclusion_time_wp, "WPExclusionTime", config);
    setConfigValue(m_leading_entries_ratio, "WPLeadingEntriesRatio", config);
    int nbins = 0;
    double xmin = 0.0, xmax = 0.0;
    if (
        !setConfigValue(nbins, "Histogram_Nbins_CD", config) ||
        !setConfigValue(xmin, "Histogram_Xmin_CD", config) ||
        !setConfigValue(xmax, "Histogram_Xmax_CD", config)
    ) return;
    m_hist_cd = std::make_unique<TH1D>((m_name + "__Histogram_CD").c_str(), (m_name + "__Histogram_CD").c_str(), nbins, xmin, xmax);
    m_hist_cd->SetDirectory(0);
    if (
        !setConfigValue(nbins, "Histogram_Nbins_WP", config) ||
        !setConfigValue(xmin, "Histogram_Xmin_WP", config) ||
        !setConfigValue(xmax, "Histogram_Xmax_WP", config)
    ) return;
    m_hist_wp = std::make_unique<TH1D>((m_name + "__Histogram_WP").c_str(), (m_name + "__Histogram_WP").c_str(), nbins, xmin, xmax);
    m_hist_wp->SetDirectory(0);
}

ParamsType ClusterMaxChargeInitializer::getOParamsType() {
    return ParamsType::SingleAcrylic;
}

bool ClusterMaxChargeInitializer::getCdInit(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable, double qthold, double& itime, vec3& ipos, vec3& fpos) {
    m_hist_cd->Reset();
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!it->used) continue;
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        m_hist_cd->Fill(it->fht, it->q);
    }   

    int max_bin_idx = m_hist_cd->GetMaximumBin();
    int leading_edge_idx = max_bin_idx;
    while (leading_edge_idx > 1 && m_hist_cd->GetBinContent(leading_edge_idx) > 0) {
        --leading_edge_idx;
    }
    double leading_edge_time = m_hist_cd->GetBinCenter(leading_edge_idx);
    double trailing_edge_time = m_hist_cd->GetBinCenter(max_bin_idx) + m_offset_final_cd;

    double leading_totq = 0.0;
    double trailing_totq = 0.0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!it->used) continue;
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (it->q < qthold) continue;
        if (leading_edge_time + m_window_min <= it-> fht && it->fht <= leading_edge_time + m_window_max) {
            leading_totq += it->q;
            ipos += (it->pos * it->q);
        }
        if (trailing_edge_time + m_window_min <= it-> fht && it->fht <= trailing_edge_time + m_window_max) {
            trailing_totq += it->q;
            fpos += (it->pos * it->q);
        }
    }
    if (leading_totq == 0.0 || trailing_totq == 0.0) {
        LogInfo << "Failed to calculate initial or final point in the CD initialization\n";
        return false;
    }
    itime = leading_edge_time;
    ipos /= leading_totq;
    fpos /= trailing_totq;
    return true;
}

bool ClusterMaxChargeInitializer::getWpInit(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable, double qthold, vec3& ipos, vec3& fpos) {
    m_hist_wp->Reset();
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!it->used) continue;
        if (!hasPmtType(*it, RecPmtType::PMT_WP)) continue;
        m_hist_wp->Fill(it->fht, it->q);
    }   

    int leading_edge_idx = m_hist_wp->GetMaximumBin();
    double leading_edge_time = m_hist_wp->GetBinCenter(leading_edge_idx);
    double leading_edge_entries = m_hist_wp->GetBinContent(leading_edge_idx);
    for (int bin = 1; bin <= m_hist_wp->GetNbinsX(); ++bin) {
        if (m_hist_wp->GetBinCenter(bin) >= leading_edge_time - m_exclusion_time_wp) break;
        if (m_hist_wp->GetBinContent(bin) >= leading_edge_entries * m_leading_entries_ratio) {
            leading_edge_idx = bin;
            leading_edge_time = m_hist_wp->GetBinCenter(leading_edge_idx);
            break;
        }
    }

    double begin_search_time = leading_edge_time + m_exclusion_time_wp;
    double trailing_edge_entries = 0.0;
    double trailing_edge_time = begin_search_time;
    for (int bin = leading_edge_idx; bin <= m_hist_wp->GetNbinsX(); ++bin) {
        if (m_hist_wp->GetBinCenter(bin) < begin_search_time) continue;
        if (m_hist_wp->GetBinContent(bin) <= trailing_edge_entries) continue;
        trailing_edge_entries = m_hist_wp->GetBinContent(bin);
        trailing_edge_time = m_hist_wp->GetBinCenter(bin);
    }

    double leading_totq = 0.0;
    double trailing_totq = 0.0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!it->used) continue;
        if (!hasPmtType(*it, RecPmtType::PMT_WP)) continue;
        if (it->q < qthold) continue;
        if (leading_edge_time + m_window_min <= it-> fht && it->fht <= leading_edge_time + m_window_max) {
            leading_totq += it->q;
            ipos += (it->pos * it->q);
        }
        if (trailing_edge_time + m_window_min <= it-> fht && it->fht <= trailing_edge_time + m_window_max) {
            trailing_totq += it->q;
            fpos += (it->pos * it->q);
        }
    }
    if (leading_totq == 0.0 || trailing_totq == 0.0) {
        LogInfo << "Failed to calculate initial or final point in the WP initialization\n";
        return false;
    }
    ipos /= leading_totq;
    fpos /= trailing_totq;
    return true;
}

bool ClusterMaxChargeInitializer::initiate(const RecPmtTable& table) {
    unsigned int count = std::count_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_20INCH) || hasPmtType(pmt, RecPmtType::PMT_WP); });
    LogDebug << count << " PMTs are used for the initialization\n";

    RecPmtTable::const_iterator it_cd = table.begin();
    RecPmtTable::const_iterator it_wp = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_CD); }).base();
    RecPmtTable::const_iterator end_wp = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_WP); }).base();
    bool has_cd = (std::distance(it_cd, it_wp) > 0);
    bool has_wp = (std::distance(it_wp, end_wp) > 0);

    bool cd_used = false;
    double itime = 0.0;
    vec3 ipos_cd, fpos_cd;
    if ( (m_mode & Mode::CD_ONLY) == Mode::CD_ONLY && has_cd ) {
        cd_used = getCdInit(it_cd, it_wp, m_qthold, itime, ipos_cd, fpos_cd);
    }
    bool wp_used = false;
    vec3 ipos_wp, fpos_wp;
    if ( (m_mode & Mode::WP_ONLY) == Mode::WP_ONLY && has_wp ) {
        wp_used = getWpInit(it_wp, end_wp, m_qthold, ipos_wp, fpos_wp);
    }

    if (!cd_used && !wp_used) {
        if ( (m_mode & Mode::CD_ONLY) == Mode::CD_ONLY && has_cd ) {
            cd_used = getCdInit(it_cd, it_wp, 0.0, itime, ipos_cd, fpos_cd);
        }
        if ( (m_mode & Mode::WP_ONLY) == Mode::WP_ONLY && has_wp ) {
            wp_used = getWpInit(it_wp, end_wp, 0.0, ipos_wp, fpos_wp);
        }
        if (!cd_used && !wp_used) {
            LogWarn << "Failed to get initialization from both CD and WP\n";
            return false;
        }
    }

    vec3 ipos = (ipos_cd * static_cast<double>(cd_used) + ipos_wp * static_cast<double>(wp_used)) / (static_cast<double>(cd_used) + static_cast<double>(wp_used));
    vec3 fpos = (fpos_cd * static_cast<double>(cd_used) + fpos_wp * static_cast<double>(wp_used)) / (static_cast<double>(cd_used) + static_cast<double>(wp_used));
    vec3 dir = unit(fpos - ipos);

    double b_half = dot(dir, ipos);
    double c = dot(ipos, ipos) - 17700.0 * 17700.0;
    double discriminant = b_half * b_half - c;
    if (discriminant >= 0.0) {
        double sqrt_discriminant = std::sqrt(discriminant);
        double t1 = -b_half - sqrt_discriminant;
        // double t2 = -b_half + sqrt_discriminant;
        ipos = ipos + t1 * dir;
    }

    m_params = std::vector<double>{itime, theta(ipos), phi(ipos), theta(dir), phi(dir)};
    return true;
}

ClusterMaxChargeInitializer::Mode operator&(const ClusterMaxChargeInitializer::Mode& a, const ClusterMaxChargeInitializer::Mode& b) {
    return static_cast<ClusterMaxChargeInitializer::Mode>(static_cast<int>(a) & static_cast<int>(b));
}