#include "initializer/fht/WaterPhaseInitializer.hpp"

#include "SniperKernel/SniperLog.h"

#include "utils/constants.hpp"

WaterPhaseInitializer::WaterPhaseInitializer(const std::string& name, double dt_f2itime, double dt, double q_ratio) : 
    Initializer<FhtMethodTag>{name},
    m_dt_f2itime{dt_f2itime},
    m_dt{dt},
    m_q_ratio{q_ratio},
    m_hist{std::make_unique<TH1D>((m_name + "__Histogram").c_str(), (m_name + "__Histogram").c_str(), 200, 0.0, 1000.0)}
{
    m_hist->SetDirectory(0);
}

void WaterPhaseInitializer::configure(const SniperJSON& config) {
    if (!config.valid()) return;
    setConfigValue(m_dt_f2itime, "MaxTimeFromExit", config);
    setConfigValue(m_dt, "EntryTimeWindow", config);
    setConfigValue(m_q_ratio, "ChargeRatio", config);
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

ParamsType WaterPhaseInitializer::getOParamsType() {
    return ParamsType::SingleCd;
}

bool WaterPhaseInitializer::isHighlyClipping(const RecPmtTable& table) {
    std::unique_ptr<TH1D> htmp = std::make_unique<TH1D>("htmp", "htmp", 100, 0.0, 1000.0);
    htmp->SetDirectory(0);
    for (const RecPmtProp& pmt : table) {
        if (!pmt.used || !hasPmtType(pmt, RecPmtType::PMT_20INCH)) continue;
        if (pmt.q < 30.0) continue;
        htmp->Fill(pmt.fht);
    }
    
    int idx_min = htmp->GetMaximumBin();
    while (htmp->GetBinContent(idx_min) > 1.0 && idx_min > 1) {
        --idx_min;
    }
    double time_min = htmp->GetBinCenter(idx_min);

    int idx_max = htmp->GetMaximumBin();
    while (htmp->GetBinContent(idx_max) > 30.0 && idx_max < htmp->GetNbinsX()) {
        ++idx_max;
    }
    double time_max = htmp->GetBinCenter(idx_max);
    
    return (time_max - time_min < 100.0);
}

bool WaterPhaseInitializer::initHighlyClipping(const RecPmtTable& table) {
    std::size_t count = std::count_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_WP); });
    bool has_wp = (count > 0ul);
    if (has_wp) {
        double time_early = std::numeric_limits<double>::infinity();
        double time_late = -std::numeric_limits<double>::infinity();

        bool found_early = false;
        bool found_late = false;

        RecPmtTable::const_iterator ftable = std::find_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_WP); });
        RecPmtTable::const_iterator ltable = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_WP); }).base();

        LogDebug << "Number of WP PMTs for initialization: " << std::distance(ftable, ltable) << '\n';

        for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
            if (!it->used || !hasPmtType(*it, RecPmtType::PMT_WP)) continue;
            if (it->q < 30.0) continue;
            time_early = std::min(time_early, it->fht);
            time_late = std::max(time_late, it->fht);
            found_early = true;
            found_late = true;
        }

        if (!found_early || !found_late) {
            for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
                if (!it->used || !hasPmtType(*it, RecPmtType::PMT_WP)) continue;
                if (it->q < 10.0) continue;
                time_early = std::min(time_early, it->fht);
                time_late = std::max(time_late, it->fht);
                found_early = true;
                found_late = true;
            }
        }

        double time_mid = 0.5 * (time_early + time_late);

        vec3 pos_early;
        vec3 pos_late;
        double totq_early = 0.0, totq_late = 0.0;

        found_early = false;
        found_late = false;

        for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
            if (!it->used || !hasPmtType(*it, RecPmtType::PMT_WP)) continue;
            if (it->q < 30.0) continue;
            if (time_early <= it->fht && it->fht <= time_mid) {
                pos_early += it->pos * it->q;
                totq_early += it->q;
                found_early = true;
            }
            if (time_mid <= it->fht && it->fht <= time_late) {
                pos_late += it->pos * it->q;
                totq_late += it->q;
                found_late = true;
            }
        }

        if (!found_early || !found_late) {
            LogWarn << "Early/Late position not found, cannot continue WpGeometryTimeTransformer" << std::endl;
            return false;
        }

        pos_early /= totq_early;
        pos_late /= totq_late;

        vec3 dir = unit(pos_late - pos_early);
        m_params = std::vector<double>{time_early, theta(pos_early), phi(pos_early), theta(dir), phi(dir)};
    }
    else {
        m_it_table.clear();
        for (RecPmtTable::const_iterator it = table.begin(); it != table.end(); ++it) {
            if (!it->used || !hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
            if (it->q < 30.0) continue;
            m_it_table.push_back(it);
        }
        if (m_it_table.empty()) {
            LogWarn << "No valid PMTs found for highly clipping initialization\n";
            return false;
        }
        double itime = std::numeric_limits<double>::infinity();
        vec3 ipos;
        double ftime = -std::numeric_limits<double>::infinity();
        vec3 fpos;
        for (RecPmtTable::const_iterator it : m_it_table) {
            if (it->fht < itime) {
                itime = it->fht;
                ipos = it->pos;
            }
            if (ftime < it->fht) {
                ftime = it->fht;
                fpos = it->pos;
            }
        }
        vec3 dir = unit(fpos - ipos);
        m_params = std::vector<double>{itime, theta(ipos), phi(ipos), theta(dir), phi(dir)};
    }
    return true;
}

bool WaterPhaseInitializer::getTable(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) {
    m_it_table.clear();
    m_hist->Reset();
    
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!it->used || !hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if ( (it->type == RecPmtType::PMT_20INCH_NNVT) && (it->q < 20.0) ) continue;
        else if ( (it->type == RecPmtType::PMT_20INCH_HAMAMATSU) && (it->q < 5.0) ) continue;
        m_hist->Fill(it->fht);
    }

    int idx = m_hist->GetMaximumBin();
    double nb_thold = 0.2 * m_hist->GetBinContent(idx);
    while (nb_thold < m_hist->GetBinContent(idx) && idx < m_hist->GetNbinsX()) {
        ++idx;
    }
    double ftime = m_hist->GetBinCenter(idx);

    double itime = std::numeric_limits<double>::infinity();
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!it->used || !hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (it->q < 30.0) continue;
        if (it->fht < ftime - m_dt_f2itime) continue;
        itime = std::min(itime, it->fht);
    }
    LogDebug << "itime: " << itime << ", ftime: " << ftime << ", ftime - f2itime: " << ftime - m_dt_f2itime << '\n';

    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!it->used || !hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (it->fht < itime || itime + m_dt < it->fht) continue;
        m_it_table.push_back(it);
    }

    if (m_it_table.empty()) {
        LogWarn << "No valid PMTs found in the time window\n";
        return false;
    }
    return true;
}

double WaterPhaseInitializer::getITime() {
    double itime = 0.0;
    double totq = 0.0;
    for (RecPmtTable::const_iterator it : m_it_table) {
        itime += it->fht * it->q;
        totq += it->q;
        
    }
    return itime / totq;
}

vec3 WaterPhaseInitializer::getIPos() {
    vec3 ipos;
    double totq = 0;
    for (RecPmtTable::const_iterator it : m_it_table) {
        ipos += it->pos * it->q;
        totq += it->q;

    }
    return ipos / totq;
}

vec3 WaterPhaseInitializer::getFPos(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable, double itime) {
    vec3 fpos;
    double totq = 0.0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!it->used || !hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (it->fht < itime + m_dt) continue;
        if (it->q < 10.0) continue;
        fpos += it->pos * it->q;
        totq += it->q;

    }
    return fpos / totq;
}

bool WaterPhaseInitializer::initiate(const RecPmtTable& table) {
    std::size_t count = std::count_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_20INCH); });
    LogDebug << count << " PMTs are used for the initialization\n";

    RecPmtTable::const_iterator ftable = std::find_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_20INCH); });
    RecPmtTable::const_iterator ltable = std::find_if(table.rbegin(), table.rend(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_20INCH); }).base();

    if (isHighlyClipping(table)) return initHighlyClipping(table);

    if (!getTable(ftable, ltable)) return false;

    double itime = getITime();
    vec3 ipos = getIPos();
    vec3 fpos = getFPos(ftable, ltable, itime);

    double ipos_theta = theta(ipos);
    double ipos_phi = phi(ipos);
    ipos = from_spherical(constants::r_cd, ipos_theta, ipos_phi);
    double fpos_theta = theta(fpos);
    double fpos_phi = phi(fpos);
    fpos = from_spherical(constants::r_cd, fpos_theta, fpos_phi);

    vec3 dir = unit(fpos - ipos);

    m_params = std::vector<double>{itime, theta(ipos), phi(ipos), theta(dir), phi(dir)};

    return true;
}