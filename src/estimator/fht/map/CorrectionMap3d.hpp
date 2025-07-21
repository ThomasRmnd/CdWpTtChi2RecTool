#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP3D_H_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP3D_H_

#include "estimator/fht/map/CorrectionMap.hpp"

#include <algorithm>

#include <TProfile3D.h>

#include "estimator/fht/map/CorrParam.hpp"

template<ParamsType _Pt>
class CorrectionMap3d : public CorrectionMap<_Pt> {

public:

    CorrectionMap3d(const std::string& name) :
        CorrectionMap<_Pt>(name)
    {}

    CorrectionMap3d(const std::string& name, const RecPmtType& pmt_type, const std::string& filename, const std::string& mapname, const std::shared_ptr<CorrParam>& param_x, const std::shared_ptr<CorrParam>& param_y, const std::shared_ptr<CorrParam>& param_z) :
        CorrectionMap<_Pt>(name, pmt_type, filename, mapname),
        m_param_x(param_x),
        m_param_y(param_y),
        m_param_z(param_z)
    {}

    ~CorrectionMap3d() override = default;

    bool initialize() override {
        if (!CorrectionMap<_Pt>::initialize()) return false;
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

    void correct(RecPmtTable& table) override {
        m_param_x->setTrack(this->m_orig, this->m_dir);
        m_param_y->setTrack(this->m_orig, this->m_dir);
        m_param_z->setTrack(this->m_orig, this->m_dir);
        for (RecPmtProp& pmt : table) {
            if (!this->checkPmtType(pmt)) continue;
            pmt.fht -= correction(pmt);
        }
        return;
    }

private:

    std::shared_ptr<CorrParam> m_param_x;
    std::shared_ptr<CorrParam> m_param_y;
    std::shared_ptr<CorrParam> m_param_z;
    TProfile3D* m_prof3d;

    int m_nbins_x, m_nbins_y, m_nbins_z;
    double m_fbin_ctr_x, m_fbin_ctr_y, m_fbin_ctr_z;
    double m_lbin_ctr_x, m_lbin_ctr_y, m_lbin_ctr_z;

    bool openCorrProfile() override {
        m_prof3d = this->m_file->template get<TProfile3D>(this->m_mapname.c_str());
        if (!m_prof3d) {
            LogError << "Cannot find TProfile3D prof3d in correction file " << this->m_filename << '\n';
            return false;
        }
        LogDebug << "TProfile3D in correction file " << this->m_filename << " opened successfully\n";

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

    int getBinX(double value) {
        int bin = m_prof3d->GetXaxis()->FindBin(value);
        if (bin == 0 || bin == m_nbins_x + 1) {
            // LogWarn << "The range of the correction map is not wide enough in the x-axis: value = " << value << ", bin = " << bin << '\n';
            bin = std::clamp(bin, 1, m_nbins_x);
        }
        return bin;
    }

    int getBinY(double value) {
        int bin = m_prof3d->GetYaxis()->FindBin(value);
        if (bin == 0 || bin == m_nbins_y + 1) {
            // LogWarn << "The range of the correction map is not wide enough in the y-axis: value = " << value << ", bin = " << bin << '\n';
            bin = std::clamp(bin, 1, m_nbins_y);
        }
        return bin;
    }

    int getBinZ(double value) {
        int bin = m_prof3d->GetZaxis()->FindBin(value);
        if (bin == 0 || bin == m_nbins_z + 1) {
            // LogWarn << "The range of the correction map is not wide enough in the z-axis: value = " << value << ", bin = " << bin << '\n';
            bin = std::clamp(bin, 1, m_nbins_z);
        }
        return bin;
    }

    double correction(const RecPmtProp& pmt) override {
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

};

template<>
class CorrectionMap3d<ParamsType::DoubleAcrylic> : public CorrectionMap<ParamsType::DoubleAcrylic> {

public:

    CorrectionMap3d(const std::string& name) :
        CorrectionMap<ParamsType::DoubleAcrylic>(name)
    {}

    CorrectionMap3d(const std::string& name, const RecPmtType& pmt_type, const std::string& filename, const std::string& mapname, const std::shared_ptr<CorrParam>& param_x, const std::shared_ptr<CorrParam>& param_y, const std::shared_ptr<CorrParam>& param_z) :
        CorrectionMap<ParamsType::DoubleAcrylic>(name, pmt_type, filename, mapname),
        m_param_x(param_x),
        m_param_y(param_y),
        m_param_z(param_z)
    {}

    ~CorrectionMap3d() override = default;

    bool initialize() override {
        if (!CorrectionMap<ParamsType::DoubleAcrylic>::initialize()) return false;
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

    void setTrack(const double* params) override {
        CorrectionMap<ParamsType::DoubleAcrylic>::setTrack(params);
        m_half_length_1 = 0.5 * m_length_1;
        m_t_end_1 = m_t_0_1 + m_length_1 / constants::c; // m_length_1 * constants::inv_c
        m_p_end_1 = m_orig_1 + m_dir * m_length_1;
        m_half_length_2 = 0.5 * m_length_2;
        m_t_end_2 = m_t_0_2 + m_length_2 / constants::c; // m_length_2 * constants::inv_c
        m_p_end_2 = m_orig_2 + m_dir * m_length_2;
    }

    void correct(RecPmtTable& table) override {
        for (RecPmtProp& pmt : table) {
            if (!checkPmtType(pmt)) continue;
            pmt.fht -= correction(pmt);
        }
        return;
    }

private:

    std::shared_ptr<CorrParam> m_param_x;
    std::shared_ptr<CorrParam> m_param_y;
    std::shared_ptr<CorrParam> m_param_z;
    TProfile3D* m_prof3d;

    int m_nbins_x, m_nbins_y, m_nbins_z;
    double m_fbin_ctr_x, m_fbin_ctr_y, m_fbin_ctr_z;
    double m_lbin_ctr_x, m_lbin_ctr_y, m_lbin_ctr_z;

    double m_half_length_1, m_half_length_2;
    double m_t_end_1, m_t_end_2;
    vec3 m_p_end_1, m_p_end_2;

    double m_fht_1, m_fht_2;
    double m_dist_orig_to_pmt_perp, m_dist_1st_light;

    bool openCorrProfile() override {
        m_prof3d = m_file->get<TProfile3D>(m_mapname.c_str());
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

    int getBinX(double value) {
        int bin = m_prof3d->GetXaxis()->FindBin(value);
        if (bin == 0 || bin == m_nbins_x + 1) {
            // LogWarn << "The range of the correction map is not wide enough in the x-axis: value = " << value << ", bin = " << bin << '\n';
            bin = std::clamp(bin, 1, m_nbins_x);
        }
        return bin;
    }

    int getBinY(double value) {
        int bin = m_prof3d->GetYaxis()->FindBin(value);
        if (bin == 0 || bin == m_nbins_y + 1) {
            // LogWarn << "The range of the correction map is not wide enough in the y-axis: value = " << value << ", bin = " << bin << '\n';
            bin = std::clamp(bin, 1, m_nbins_y);
        }
        return bin;
    }

    int getBinZ(double value) {
        int bin = m_prof3d->GetZaxis()->FindBin(value);
        if (bin == 0 || bin == m_nbins_z + 1) {
            // LogWarn << "The range of the correction map is not wide enough in the z-axis: value = " << value << ", bin = " << bin << '\n';
            bin = std::clamp(bin, 1, m_nbins_z);
        }
        return bin;
    }

    double correction(const RecPmtProp& pmt) override {
        m_dist_orig_to_pmt_perp = m_half_length_1 + dot(m_dir, pmt.pos);
        m_dist_1st_light = m_dist_orig_to_pmt_perp - mag(pmt.pos - (m_orig_1 + m_dist_orig_to_pmt_perp * m_dir)) * constants::inv_tan_cherenkov_ls;
        if (m_dist_1st_light < 0.0) m_fht_1 = m_t_0_1 + constants::inv_c_ls * mag(pmt.pos - m_orig_1);
        else if (m_length_1 < m_dist_1st_light) m_fht_1 = m_t_end_1 + constants::inv_c_ls * mag(pmt.pos - m_p_end_1);
        else m_fht_1 = m_t_0_1 + constants::one_minus_n_ls_2_div_c * m_dist_1st_light + constants::n_ls_2_div_c * m_dist_orig_to_pmt_perp;

        m_dist_orig_to_pmt_perp = m_half_length_2 + dot(m_dir, pmt.pos);
        m_dist_1st_light = m_dist_orig_to_pmt_perp - mag(pmt.pos - (m_orig_2 + m_dist_orig_to_pmt_perp * m_dir)) * constants::inv_tan_cherenkov_ls;
        if (m_dist_1st_light < 0.0) m_fht_2 = m_t_0_2 + constants::inv_c_ls * mag(pmt.pos - m_orig_2);
        else if (m_length_2 < m_dist_1st_light) m_fht_2 = m_t_end_2 + constants::inv_c_ls * mag(pmt.pos - m_p_end_2);
        else m_fht_2 = m_t_0_2 + constants::one_minus_n_ls_2_div_c * m_dist_1st_light + constants::n_ls_2_div_c * m_dist_orig_to_pmt_perp;

        if (m_fht_1 < m_fht_2) {
            m_param_x->setTrack(m_orig_1, m_dir);
            m_param_y->setTrack(m_orig_1, m_dir);
            m_param_z->setTrack(m_orig_1, m_dir);
        }
        else {
            m_param_x->setTrack(m_orig_2, m_dir);
            m_param_y->setTrack(m_orig_2, m_dir);
            m_param_z->setTrack(m_orig_2, m_dir);
        }

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

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP3D_H_