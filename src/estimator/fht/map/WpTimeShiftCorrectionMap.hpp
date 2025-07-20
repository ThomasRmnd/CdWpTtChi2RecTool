#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_

#include "estimator/fht/map/CorrectionMap.hpp"

#include <memory>

#include <TProfile2D.h>

#include "utils/Geometry.hpp"

template<ParamsType _Pt>
class WpTimeShiftCorrectionMap : public CorrectionMap<_Pt> {

public:

    WpTimeShiftCorrectionMap(const std::string& name) :
        CorrectionMap<_Pt>(name)
    {};

    WpTimeShiftCorrectionMap(const std::string& name, const std::string& filename, const std::string& mapname) :
        CorrectionMap<_Pt>(name, RecPmtType::PMT_WP, filename, mapname)
    {};

    ~WpTimeShiftCorrectionMap() override = default;

    void operator()(RecPmtTable& table) override {
    
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
            if (!this->check(pmt)) continue;
            signed_pos = dot(this->m_dir, pmt.pos);
            if (signed_pos <= 0.0) continue; // we only correct 'bottom' PMTs
            pmt.fht += correction(pmt);
        }
        return;
    };

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
    };

    int getBin(TAxis* axis, double value) {
        int bin = axis->FindBin(value);
        if (bin == 0 || bin == axis->GetNbins() + 1) {
            // LogWarn << "The range of the correction map is not wide enough: value = " << value << ", bin = " << bin << std::endl;
            bin = std::clamp(bin, 1, axis->GetNbins());
        }
        return bin;
    };

    double correction(const RecPmtProp&) {
        return m_shift;
    };

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_