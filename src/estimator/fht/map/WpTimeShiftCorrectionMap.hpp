#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_

#include <memory>

#include <TProfile2D.h>

#include "estimator/fht/map/CorrectionMap.hpp"
#include "utils/Geometry.hpp"

template<typename _ParamsTag>
class WpTimeShiftCorrectionMap : public CorrectionMap, public TrackSetter<_ParamsTag> {

public:

    WpTimeShiftCorrectionMap(const std::string& name, const std::string& filename, const std::string& mapname) :
        CorrectionMap(name, RecPmtType::PMT_WP, filename, mapname)
    {}

    ~WpTimeShiftCorrectionMap() override = default;

    ParamsType getIParamsType() override {
        return ParamsTraits<_ParamsTag>::type;
    }

    void setTrack(const double* params) override {
        TrackSetter<_ParamsTag>::setTrack(params);
        auto [d_i_wp, d_o_wp] = computeWpHits(this->m_orig, this->m_dir);
        m_p_end_wp = this->m_orig + d_o_wp * this->m_dir;
    }

    void correct(RecPmtTable::iterator first, RecPmtTable::iterator last, const double* params) override {
        this->setTrack(params);

        double x = theta(m_p_end_wp);
        double y = theta(this->m_dir);

        int bin_x = this->getBin(m_prof->GetXaxis(), x);
        int bin_y = this->getBin(m_prof->GetYaxis(), y);

        m_shift = m_prof->GetBinContent(bin_x, bin_y);

        LogDebug << "x = " << x << " y = " << y << std::endl;
        LogDebug << "WP Time Shift: " << m_shift << std::endl;

        double signed_pos = 0.0;
        for (RecPmtTable::iterator it = first; it != last; ++it) {
            if (!checkPmtType(*it)) continue;
            signed_pos = dot(this->m_dir, it->pos);
            if (signed_pos <= 0.0) continue; // we only correct 'bottom' PMTs
            it->fht += correction(*it);
        }
        return;
    }

private:

    TProfile2D* m_prof;
    double m_shift;

    vec3 m_p_end_wp;

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
class WpTimeShiftCorrectionMap<DoubleAcrylicParamsTag> : public CorrectionMap, public TrackSetter<DoubleAcrylicParamsTag> {

public:

    WpTimeShiftCorrectionMap(const std::string& name, const std::string& filename, const std::string& mapname) :
        CorrectionMap(name, RecPmtType::PMT_WP, filename, mapname)
    {}

    ~WpTimeShiftCorrectionMap() override = default;

    ParamsType getIParamsType() override {
        return ParamsTraits<DoubleAcrylicParamsTag>::type;
    }

    void setTrack(const double* params) override {
        TrackSetter<DoubleAcrylicParamsTag>::setTrack(params);
        std::tie(m_d_i_wp_1, m_d_o_wp_1) = computeWpHits(m_orig_1, m_dir);
        m_p_end_wp_1 = m_orig_1 + m_d_o_wp_1 * m_dir;
        std::tie(m_d_i_wp_2, m_d_o_wp_2) = computeWpHits(m_orig_2, m_dir);
        m_p_end_wp_2 = m_orig_2 + m_d_o_wp_2 * m_dir;
    }

    void correct(RecPmtTable::iterator first, RecPmtTable::iterator last, const double* params) override {
        setTrack(params);
        double signed_pos = 0.0;
        for (RecPmtTable::iterator it = first; it != last; ++it) {
            if (!checkPmtType(*it)) continue;
            signed_pos = dot(m_dir, it->pos);
            if (signed_pos <= 0.0) continue; // we only correct 'bottom' PMTs
            it->fht += correction(*it);
        }
        return;
    }

private:

    TProfile2D* m_prof;

    double m_d_i_wp_1, m_d_o_wp_1;
    vec3 m_p_end_wp_1;
    double m_d_i_wp_2, m_d_o_wp_2;
    vec3 m_p_end_wp_2;

    bool openCorrProfile() override {
        m_prof = m_file->template get<TProfile2D>(m_mapname.c_str());
        if (!m_prof) {
            LogError << "Cannot find TProfile3D prof3d in correction file " << m_filename << std::endl;
            return false;
        }
        LogDebug << "WP Profile in correction file " << m_filename << " opened successfully." << std::endl;
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
        double fht_1 = m_t_0_1 + m_d_o_wp_1 * constants::inv_c + mag(pmt.pos - m_p_end_wp_1) * constants::inv_c_w - 10.0; // 10 is a constant, maybe because of the time emission
        double fht_2 = m_t_0_2 + m_d_o_wp_2 * constants::inv_c + mag(pmt.pos - m_p_end_wp_2) * constants::inv_c_w - 10.0;

        vec3 end_wp = (fht_1 < fht_2 ? m_p_end_wp_1 : m_p_end_wp_2);

        double x = theta(end_wp);
        double y = theta(m_dir);

        int bin_x = getBin(m_prof->GetXaxis(), x);
        int bin_y = getBin(m_prof->GetYaxis(), y);

        return m_prof->GetBinContent(bin_x, bin_y);
    }

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_WPTIMESHIFTCORRECTIONMAP_HPP_