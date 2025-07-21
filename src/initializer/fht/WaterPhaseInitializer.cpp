#include "initializer/fht/WaterPhaseInitializer.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

#include "utils/constants.hpp"

DECLARE_TOOL(WaterPhaseInitializer);

WaterPhaseInitializer::WaterPhaseInitializer(const std::string& name) : 
    Initializer<FhtMethodTag>{name},
    m_dt_f2itime{200.0},
    m_dt{20.0},
    m_q_ratio{0.5},
    m_hist{std::make_unique<TH1D>("h__WaterPhaseInitializer", "h__WaterPhaseInitializer", 200, 0.0, 1000.0)}
{
    m_hist->SetDirectory(0);
}

WaterPhaseInitializer::WaterPhaseInitializer(const std::string& name, double dt_f2itime, double dt, double q_ratio) : 
    Initializer<FhtMethodTag>{name},
    m_dt_f2itime{dt_f2itime},
    m_dt{dt},
    m_q_ratio{q_ratio},
    m_hist{std::make_unique<TH1D>("hist", "hist", 200, 0.0, 1000.0)}
{
    m_hist->SetDirectory(0);
}

bool WaterPhaseInitializer::getTable(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) {
    m_it_table.clear();
    m_hist->Reset();
    
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
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
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (it->q < 30.0) continue;
        if (it->fht < ftime - m_dt_f2itime) continue;
        itime = std::min(itime, it->fht);
    }
    LogDebug << "itime: " << itime << ", ftime: " << ftime << ", ftime - f2itime: " << ftime - m_dt_f2itime << '\n';

    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
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
    double qmax = 0.0;
    vec3 posmax;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (it->fht < itime + m_dt) continue;
        if (qmax < it->q) {
            qmax = it->q;
            posmax = it->pos;
        }

    }
    LogDebug << "qmax = " << qmax << '\n';

    vec3 fpos;
    double totq = 0.0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (!hasPmtType(*it, RecPmtType::PMT_20INCH)) continue;
        if (it->fht < itime + m_dt) continue;
        if (it->q < qmax * m_q_ratio) continue;
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

    if (!getTable(ftable, ltable)) return false;

    double itime = getITime();
    vec3 ipos = getIPos();
    vec3 fpos = getFPos(ftable, ltable, itime);

#define __CDWPTTCHI2RECTOOL_WATERPHASEINITIALIZER_SIMPLE_OUTPUT_PARAMS__

#ifdef __CDWPTTCHI2RECTOOL_WATERPHASEINITIALIZER_SIMPLE_OUTPUT_PARAMS__
    double ipos_theta = theta(ipos);
    double ipos_phi = phi(ipos);
    ipos = from_spherical(constants::r_cd, ipos_theta, ipos_phi);
    double fpos_theta = theta(fpos);
    double fpos_phi = phi(fpos);
    fpos = from_spherical(constants::r_cd, fpos_theta, fpos_phi);
#endif // __CDWPTTCHI2RECTOOL_WATERPHASEINITIALIZER_SIMPLE_OUTPUT_PARAMS__

    vec3 dir = unit(fpos - ipos);

#ifndef __CDWPTTCHI2RECTOOL_WATERPHASEINITIALIZER_SIMPLE_OUTPUT_PARAMS__
    double _b = dot(dir, ipos);
    double _c = mag2(ipos) - constants::r_cd * constants::r_cd;

    double _delta = _b * _b - _c;
    if (_delta < 0.0) {
        ipos_theta = theta(ipos);
        ipos_phi = phi(ipos);
        ipos = from_spherical(constants::r_cd, ipos_theta, ipos_phi);
        fpos_theta = theta(fpos);
        fpos_phi = phi(fpos);
        fpos = from_spherical(constants::r_cd, fpos_theta, fpos_phi);
    }
    else {
        double _d2 = -_b + std::sqrt(_delta);
        fpos = ipos + dir * _d2;
        double _d1 = -_b - std::sqrt(_delta);
        ipos = ipos + dir * _d1;
    }

    dir = unit(fpos - ipos);
#endif // __CDWPTTCHI2RECTOOL_WATERPHASEINITIALIZER_SIMPLE_OUTPUT_PARAMS__

    m_params = std::vector<double>{itime, theta(ipos), phi(ipos), theta(dir), phi(dir)};

    return true;
}

ParamsType WaterPhaseInitializer::getOParamsType() {
    return ParamsType::SingleCd;
}