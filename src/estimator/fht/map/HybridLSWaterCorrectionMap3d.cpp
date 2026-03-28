#include "estimator/fht/map/HybridLSWaterCorrectionMap3d.hpp"

#include <algorithm>
#include <limits>

#include "utils/constants.hpp"

HybridLSWaterCorrectionMap3d::HybridLSWaterCorrectionMap3d(const std::string& name, const RecPmtType& pmt_type, const std::string& filename, const std::string& mapname, const std::shared_ptr<CorrParam>& param_x, const std::shared_ptr<CorrParam>& param_y, const std::shared_ptr<CorrParam>& param_z) :
    CorrectionMap{name, pmt_type, filename, mapname},
    m_param_x{param_x},
    m_param_y{param_y},
    m_param_z{param_z}
{}

ParamsType HybridLSWaterCorrectionMap3d::getIParamsType() {
    return ParamsType::SingleCd;
}

bool HybridLSWaterCorrectionMap3d::initialize() {
    if (!CorrectionMap::initialize()) return false;
    if (!m_param_x) {
        LogError << "CorrParam x is not set\n";
        return false;
    }
    if (!m_param_y) {
        LogError << "CorrParam y is not set\n";
        return false;
    }
    if (!m_param_z) {
        LogError << "CorrParam z is not set\n";
        return false;
    }
    return true;
}

void HybridLSWaterCorrectionMap3d::correct(RecPmtTable::iterator first, RecPmtTable::iterator last, const double* params) {
    setTrack(params);
    double b_half = dot(m_dir, m_orig);
    double c = mag2(m_orig) - constants::r_acrylic * constants::r_acrylic;
    double disc_fourth = b_half * b_half - c;
    if (disc_fourth < 0.0) return;
    double d = -b_half - std::sqrt(disc_fourth);
    vec3 orig_ = m_orig + d * m_dir;
    const double params_[5] = {params[0] + d / constants::c, theta(orig_), phi(orig_), params[3], params[4]};
    setTrack(params_);

    m_param_x->setTrack(m_orig, m_dir);
    m_param_y->setTrack(m_orig, m_dir);
    m_param_z->setTrack(m_orig, m_dir);
    for (RecPmtTable::iterator it = first; it != last; ++it) {
        if (!checkPmtType(*it)) continue;
        it->fht -= correction(*it);
    }
    return;
}

bool HybridLSWaterCorrectionMap3d::openCorrProfile() {
    m_prof3d = m_file->template get<TProfile3D>(m_mapname.c_str());
    if (!m_prof3d) {
        LogError << "Cannot find TProfile3D prof3d in correction file " << m_filename << '\n';
        return false;
    }
    LogDebug << "TProfile3D in correction file " << m_filename << " opened successfully\n";

    m_nbins_x = m_prof3d->GetXaxis()->GetNbins();
    m_nbins_y = m_prof3d->GetYaxis()->GetNbins();
    m_nbins_z = m_prof3d->GetZaxis()->GetNbins();
    m_fbin_ctr_x = m_prof3d->GetXaxis()->GetBinCenter(1);
    m_fbin_ctr_y = m_prof3d->GetYaxis()->GetBinCenter(1);
    m_fbin_ctr_z = m_prof3d->GetZaxis()->GetBinCenter(1);
    m_lbin_ctr_x = m_prof3d->GetXaxis()->GetBinCenter(m_prof3d->GetXaxis()->GetNbins());
    m_lbin_ctr_y = m_prof3d->GetYaxis()->GetBinCenter(m_prof3d->GetYaxis()->GetNbins());
    m_lbin_ctr_z = m_prof3d->GetZaxis()->GetBinCenter(m_prof3d->GetZaxis()->GetNbins());
    
    return true;
}

int HybridLSWaterCorrectionMap3d::getBinX(double value) {
    int bin = m_prof3d->GetXaxis()->FindBin(value);
    if (bin == 0 || bin == m_nbins_x + 1) {
        // LogWarn << "The range of the correction map is not wide enough in the x-axis: value = " << value << ", bin = " << bin << '\n';
        bin = std::clamp(bin, 1, m_nbins_x);
    }
    return bin;
}

int HybridLSWaterCorrectionMap3d::getBinY(double value) {
    int bin = m_prof3d->GetYaxis()->FindBin(value);
    if (bin == 0 || bin == m_nbins_y + 1) {
        // LogWarn << "The range of the correction map is not wide enough in the y-axis: value = " << value << ", bin = " << bin << '\n';
        bin = std::clamp(bin, 1, m_nbins_y);
    }
    return bin;
}

int HybridLSWaterCorrectionMap3d::getBinZ(double value) {
    int bin = m_prof3d->GetZaxis()->FindBin(value);
    if (bin == 0 || bin == m_nbins_z + 1) {
        // LogWarn << "The range of the correction map is not wide enough in the z-axis: value = " << value << ", bin = " << bin << '\n';
        bin = std::clamp(bin, 1, m_nbins_z);
    }
    return bin;
}

double HybridLSWaterCorrectionMap3d::correction(const RecPmtProp& pmt) {
    double corr = 0.0;

    double x = m_param_x->calculate(pmt);
    double y = m_param_y->calculate(pmt);
    double z = m_param_z->calculate(pmt);
    
    if ( (m_fbin_ctr_x < x && x < m_lbin_ctr_x) && (m_fbin_ctr_y < y && y < m_lbin_ctr_y) && (m_fbin_ctr_z < z && z < m_lbin_ctr_z) ) {
        corr = m_prof3d->Interpolate(x, y, z);
    }
    else { // Out of range, cannot interpolate
        int x_bin = getBinX(x);
        int y_bin = getBinY(y);
        int z_bin = getBinZ(z);

        corr = m_prof3d->GetBinContent(x_bin, y_bin, z_bin);
    }
    
    return corr;
}