#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_

#include "estimator/fht/map/CorrectionMap.hpp"

#include <memory>

#include <TProfile2D.h>

#include "utils/Geometry.hpp"

template<typename _ParamsTag>
class WpTimeShiftCorrectionMap : public CorrectionMap<_ParamsTag> {

public:

    WpTimeShiftCorrectionMap(const std::string& name, const std::string& filename, const std::string& mapname) :
        CorrectionMap<_ParamsTag>(name, RecPmtType::PMT_WP, filename, mapname)
    {}

    ~WpTimeShiftCorrectionMap() override = default;

    void correct(RecPmtTable& table) override {
    
        double d_i_wp = 0.0, d_o_wp = 0.0;
        if (!trajectoryHitCylinder(21750.0, 44000.0, this->m_orig, this->m_dir, d_i_wp, d_o_wp)) {
            LogError << "Cannot project points into the WP." << std::endl;
            LogError << "Orig: " << this->m_orig.x << " " << this->m_orig.y << " " << this->m_orig.z << std::endl;
            LogError << "Dir: " << this->m_dir.x << " " << this->m_dir.y << " " << this->m_dir.z << std::endl;
            return;
        }

        vec3 end_wp = this->m_orig + d_o_wp * this->m_dir;

        double x = theta(end_wp);
        double y = theta(this->m_dir);

        int bin_x = this->getBin(m_prof->GetXaxis(), x);
        int bin_y = this->getBin(m_prof->GetYaxis(), y);

        m_shift = m_prof->GetBinContent(bin_x, bin_y);

        LogDebug << "x = " << x << " y = " << y << std::endl;
        LogDebug << "WP Time Shift: " << m_shift << std::endl;

        double signed_pos = 0.0;
        for (RecPmtProp& pmt : table) {
            if (!this->checkPmtType(pmt)) continue;
            signed_pos = dot(this->m_dir, pmt.pos);
            if (signed_pos <= 0.0) continue; // we only correct 'bottom' PMTs
            pmt.fht += correction(pmt);
        }
        return;
    }

private:

    TProfile2D* m_prof;
    double m_shift;

    bool openCorrProfile() override {
        m_prof = this->m_file->template get<TProfile2D>(this->m_mapname.c_str());
        if (!m_prof) {
            LogError << "Cannot find TProfile3D prof3d in correction file " << this->m_filename << std::endl;
            return false;
        }
        LogDebug << "WP Profile in correction file " << this->m_filename << " opened successfully." << std::endl;
        return true;
    }

    int getBin(TAxis* axis, double value) {
        int bin = axis->FindBin(value);
        if (bin == 0 || bin == axis->GetNbins() + 1) {
            // LogWarn << "The range of the correction map is not wide enough: value = " << value << ", bin = " << bin << std::endl;
            bin = std::clamp(bin, 1, axis->GetNbins());
        }
        return bin;
    }

    double correction(const RecPmtProp&) {
        return m_shift;
    }

};

template<>
class WpTimeShiftCorrectionMap<DoubleAcrylicParamsTag> : public CorrectionMap<DoubleAcrylicParamsTag> {

public:

    WpTimeShiftCorrectionMap(const std::string& name, const std::string& filename, const std::string& mapname) :
        CorrectionMap<DoubleAcrylicParamsTag>(name, RecPmtType::PMT_WP, filename, mapname)
    {}

    ~WpTimeShiftCorrectionMap() override = default;

    void correct(RecPmtTable& table) override {
        double signed_pos = 0.0;
        for (RecPmtProp& pmt : table) {
            if (!checkPmtType(pmt)) continue;
            signed_pos = dot(m_dir, pmt.pos);
            if (signed_pos <= 0.0) continue; // we only correct 'bottom' PMTs
            pmt.fht += correction(pmt);
        }
        return;
    }

private:

    TProfile2D* m_prof;

    bool openCorrProfile() override {
        m_prof = this->m_file->template get<TProfile2D>(this->m_mapname.c_str());
        if (!m_prof) {
            LogError << "Cannot find TProfile3D prof3d in correction file " << this->m_filename << std::endl;
            return false;
        }
        LogDebug << "WP Profile in correction file " << this->m_filename << " opened successfully." << std::endl;
        return true;
    }

    int getBin(TAxis* axis, double value) {
        int bin = axis->FindBin(value);
        if (bin == 0 || bin == axis->GetNbins() + 1) {
            // LogWarn << "The range of the correction map is not wide enough: value = " << value << ", bin = " << bin << std::endl;
            bin = std::clamp(bin, 1, axis->GetNbins());
        }
        return bin;
    }

    double correction(const RecPmtProp& pmt) {
        double d_i_wp_1 = 0.0, d_o_wp_1 = 0.0;
        double d_i_wp_2 = 0.0, d_o_wp_2 = 0.0;
        if (
            !trajectoryHitCylinder(21750.0, 44000.0, m_orig_1, m_dir, d_i_wp_1, d_o_wp_1) ||
            !trajectoryHitCylinder(21750.0, 44000.0, m_orig_2, m_dir, d_i_wp_2, d_o_wp_2)
        ) { // need to be changed for a analytical solution
            LogError << "Cannot project points into the WP\n";
            LogError << "Orig: " << m_orig_1.x << ' ' << m_orig_1.y << ' ' << m_orig_1.z << '\n';
            LogError << "Orig: " << m_orig_2.x << ' ' << m_orig_2.y << ' ' << m_orig_2.z << '\n';
            LogError << "Dir: " << m_dir.x << ' ' << m_dir.y << ' ' << m_dir.z << '\n';
            std::cin.get();
            return 0.0;
        }

        vec3 end_wp_1 = m_orig_1 + m_dir * d_o_wp_1;
        vec3 end_wp_2 = m_orig_2 + m_dir * d_o_wp_2;

        double signed_pos = dot(m_dir, pmt.pos);

        double fht_1 = 0.0, fht_2 = 0.0;

        vec3 p_pmt_proj = m_orig_1 + dot(m_dir, pmt.pos - m_orig_1) * m_dir;
        double dt_pmt_proj = dot(p_pmt_proj - m_orig_1, m_dir) * constants::inv_c;
        double dt_1st_light = dt_pmt_proj - mag(pmt.pos - p_pmt_proj) * constants::inv_c_tan_cherenkov_w;
        double d_1st_light = constants::c * dt_1st_light;
        vec3 p_1st_light = m_orig_1 + d_1st_light * this->m_dir;
    
        fht_1 = m_t_0_1 + dt_1st_light + mag(pmt.pos - p_1st_light) * constants::inv_c_w;

        p_pmt_proj = m_orig_2 + dot(m_dir, pmt.pos - m_orig_2) * m_dir;
        dt_pmt_proj = dot(p_pmt_proj - m_orig_2, m_dir) * constants::inv_c;
        dt_1st_light = dt_pmt_proj - mag(pmt.pos - p_pmt_proj) * constants::inv_c_tan_cherenkov_w;
        d_1st_light = constants::c * dt_1st_light;
        p_1st_light = m_orig_2 + d_1st_light * this->m_dir;

        fht_2 = m_t_0_2 + dt_1st_light + mag(pmt.pos - p_1st_light) * constants::inv_c_w;

        if (0 < signed_pos) {
            fht_1 = m_t_0_1 + d_o_wp_1 * constants::inv_c + mag(pmt.pos - end_wp_1) * constants::inv_c_w /* -10.0 */; // 10 is a constant, maybe because of the time emission
            fht_2 = m_t_0_2 + d_o_wp_2 * constants::inv_c + mag(pmt.pos - end_wp_2) * constants::inv_c_w /* -10.0 */;
        }

        vec3 end_wp = (fht_1 < fht_2 ? end_wp_1 : end_wp_2);

        double x = theta(end_wp);
        double y = theta(m_dir);

        int bin_x = getBin(m_prof->GetXaxis(), x);
        int bin_y = getBin(m_prof->GetYaxis(), y);

        return m_prof->GetBinContent(bin_x, bin_y);
    }

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_