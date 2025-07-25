#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_WPFHT_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_WPFHT_HPP_

#include "SniperKernel/SniperLog.h"

#include "utils/constants.hpp"
#include "utils/Geometry.hpp"
#include "utils/RecPmtProp.hpp"
#include "utils/TrackParams.hpp"

/**
 * @class WpFht
 * @tparam _ParamsTag Method tag
 * 
 * @brief Base class to calculate expected FHT for WP PMT
 */
template<typename _ParamsType>
class WpFht : public TrackSetter<_ParamsType> {

    static_assert(std::is_base_of<ParamsTag, _ParamsTag>::value, "Tag must derive from ParamsTag");

public:

    virtual ~WpFht() override = default;

    /**
     * @brief Calculate the expected FHT
     * 
     * @param pmt PMT information
     * 
     * @return FHT
     */
    virtual double calculate(const RecPmtProp& pmt) = 0;

};

/**
 * @class NoDiffusionWpFht
 * @tparam _ParamsTag Method tag
 * 
 * @brief Derived class to calculate expected FHT for WP PMT
 * This calculation do not take into account diffusions produced by the Tyvek
 */
template<typename _ParamsType>
class NoDiffusionWpFht final : public WpFht<_ParamsType> {

public:

    ~NoDiffusionWpFht() final override = default;

    /**
     * @brief Calculate the expected FHT
     * 
     * @param pmt PMT information
     * 
     * @return FHT
     */
    double calculate(const RecPmtProp& pmt) {

        // TODO when finished checking: Following this folloing paragraph in the setTrack method {

        double d_i_wp = 0.0, d_o_wp = 0.0;
        if (!trajectoryHitCylinder(21750.0, 44000.0, this->m_orig, this->m_dir, d_i_wp, d_o_wp)) { // need to be changed for a analytical solution
            LogError << "Cannot project points into the WP\n";
            LogError << "Orig: " << this->m_orig.x << ' ' << this->m_orig.y << ' ' << this->m_orig.z << '\n';
            LogError << "Dir: " << this->m_dir.x << ' ' << this->m_dir.y << ' ' << this->m_dir.z << '\n';
            std::cin.get();
            return 0.0;
        }

        vec3 end_wp = this->m_orig + this->m_dir * d_o_wp;

        // } end paragraph

        // Cherenkov

        m_p_pmt_proj = this->m_orig + dot(this->m_dir, pmt.pos - this->m_orig) * this->m_dir;
        double dt_pmt_proj = dot(m_p_pmt_proj - this->m_orig, this->m_dir) * constants::inv_c;
        double dt_1st_light = dt_pmt_proj - mag(pmt.pos - m_p_pmt_proj) * constants::inv_c_tan_cherenkov_w;
        double d_1st_light = constants::c * dt_1st_light;
        m_p_1st_light = this->m_orig + d_1st_light * this->m_dir;
    
        double t_cherenkov = this->m_t_0 + dt_1st_light + mag(pmt.pos - m_p_1st_light) * constants::inv_c_w;

        double signed_pos = dot(this->m_dir, pmt.pos);

        if (0 < signed_pos) {
            /* if (this->m_length <= d_1st_light && d_1st_light <= d_o_wp) {
                vec3 p0 = m_p_1st_light;
                vec3 d0 = unit(pmt.pos - m_p_1st_light);
                double _b = dot(p0, d0);
                double _c = mag2(p0) - 20050.0 * 20050.0;
                double _d = _b * _b - _c;
                if (_d > 0.0) {
                    double _t = -_b - std::sqrt(_d);
                    if (_t < mag(pmt.pos - m_p_1st_light)) return t_cherenkov;
                }
                else return t_cherenkov;
            } */

            double t_back = this->m_t_0 + d_o_wp * constants::inv_c + mag(pmt.pos - end_wp) * constants::inv_c_w /* -10.0 */; // 10 is a constant, maybe because of the time emission
            return t_back;

            /* double dist_to_trk = mag(cross(m_dir, pmt.pos - m_orig));
            if (dist_to_trk <= 5000.0 && ( std::abs(pmt.fht - t_cherenkov) < std::abs(pmt.fht - t_back) )) {
                return t_cherenkov;
            }
            else {
                return t_back;
            } */
        }

        else { // signed_pos <= 0
            return t_cherenkov;
        }

    };

private:

    vec3 m_p_pmt_proj, m_p_1st_light;

    double m_dist_1st_light;
    double m_dist_orig_to_pmt_perp;

};

/**
 * @class NoDiffusionWpFht<DoubleAcrylicParamsTag>
 * 
 * @brief Derived class to calculate expected FHT for WP PMT, specialization with DoubleAcrylicParamsTag
 * This calculation do not take into account diffusions produced by the Tyvek
 */
template<>
class NoDiffusionWpFht<DoubleAcrylicParamsTag> final : public WpFht<DoubleAcrylicParamsTag> {

public:

    ~NoDiffusionWpFht() final override = default;

    /**
     * @brief Calculate the expected FHT
     * 
     * @param pmt PMT information
     * 
     * @return FHT
     */
    double calculate(const RecPmtProp& pmt) {

        // TODO when finished checking: Following this folloing paragraph in the setTrack method {

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

        // } end paragraph

        // Cherenkov

        double signed_pos = dot(m_dir, pmt.pos);

        double fht_1 = 0.0, fht_2 = 0.0;

        m_p_pmt_proj = m_orig_1 + dot(m_dir, pmt.pos - m_orig_1) * m_dir;
        double dt_pmt_proj = dot(m_p_pmt_proj - m_orig_1, m_dir) * constants::inv_c;
        double dt_1st_light = dt_pmt_proj - mag(pmt.pos - m_p_pmt_proj) * constants::inv_c_tan_cherenkov_w;
        double d_1st_light = constants::c * dt_1st_light;
        m_p_1st_light = m_orig_1 + d_1st_light * this->m_dir;
    
        fht_1 = m_t_0_1 + dt_1st_light + mag(pmt.pos - m_p_1st_light) * constants::inv_c_w;

        m_p_pmt_proj = m_orig_2 + dot(m_dir, pmt.pos - m_orig_2) * m_dir;
        dt_pmt_proj = dot(m_p_pmt_proj - m_orig_2, m_dir) * constants::inv_c;
        dt_1st_light = dt_pmt_proj - mag(pmt.pos - m_p_pmt_proj) * constants::inv_c_tan_cherenkov_w;
        d_1st_light = constants::c * dt_1st_light;
        m_p_1st_light = m_orig_2 + d_1st_light * this->m_dir;

        fht_2 = m_t_0_2 + dt_1st_light + mag(pmt.pos - m_p_1st_light) * constants::inv_c_w;

        if (0 < signed_pos) {
            fht_1 = m_t_0_1 + d_o_wp_1 * constants::inv_c + mag(pmt.pos - end_wp_1) * constants::inv_c_w /* -10.0 */; // 10 is a constant, maybe because of the time emission
            fht_2 = m_t_0_2 + d_o_wp_2 * constants::inv_c + mag(pmt.pos - end_wp_2) * constants::inv_c_w /* -10.0 */;
        }
        return std::min(fht_1, fht_2);
    };

private:

    vec3 m_p_pmt_proj, m_p_1st_light;

    double m_dist_1st_light;
    double m_dist_orig_to_pmt_perp;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_WPFHT_HPP_