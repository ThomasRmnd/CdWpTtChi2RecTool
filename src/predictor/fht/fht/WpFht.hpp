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
template<typename _ParamsTag>
class WpFht : public TrackSetter<_ParamsTag> {

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
template<typename _ParamsTag>
class NoDiffusionWpFht final : public WpFht<_ParamsTag> {

public:

    ~NoDiffusionWpFht() final override = default;

    void setTrack(const double* params) final override {
        WpFht<_ParamsTag>::setTrack(params);
        std::tie(m_d_i_wp, m_d_o_wp) = computeWpHits(this->m_orig, this->m_dir);
        m_p_end_wp = this->m_orig + this->m_dir * m_d_o_wp;
    }

    /**
     * @brief Calculate the expected FHT
     * 
     * @param pmt PMT information
     * 
     * @return FHT
     */
    double calculate(const RecPmtProp& pmt) {
        // TODO: optimize
        m_p_pmt_proj = this->m_orig + dot(this->m_dir, pmt.pos - this->m_orig) * this->m_dir;
        double dt_pmt_proj = dot(m_p_pmt_proj - this->m_orig, this->m_dir) * constants::inv_c;
        double dt_1st_light = dt_pmt_proj - mag(pmt.pos - m_p_pmt_proj) * constants::inv_c_tan_cherenkov_w;
        double d_1st_light = constants::c * dt_1st_light;
        m_p_1st_light = this->m_orig + d_1st_light * this->m_dir;
    
        double t_cherenkov = this->m_t_0 + dt_1st_light + mag(pmt.pos - m_p_1st_light) * constants::inv_c_w;

        double signed_pos = dot(this->m_dir, pmt.pos);

        if (0 < signed_pos) {
            return this->m_t_0 + m_d_o_wp * constants::inv_c + mag(pmt.pos - m_p_end_wp) * constants::inv_c_w - 10.0; // 10 is an arbitrary constant, leads to better performances ==> maybe because of the Tyvek diffusion
        }

        return t_cherenkov;
    };

private:

    double m_d_i_wp, m_d_o_wp;
    vec3 m_p_end_wp;

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

    void setTrack(const double* params) final override {
        WpFht<DoubleAcrylicParamsTag>::setTrack(params);
        std::tie(m_d_i_wp_1, m_d_o_wp_1) = computeWpHits(m_orig_1, m_dir);
        std::tie(m_d_i_wp_2, m_d_o_wp_2) = computeWpHits(m_orig_2, m_dir);
        m_p_end_wp_1 = m_orig_1 + m_dir * m_d_o_wp_1;
        m_p_end_wp_2 = m_orig_2 + m_dir * m_d_o_wp_2;
    }

    /**
     * @brief Calculate the expected FHT
     * 
     * @param pmt PMT information
     * 
     * @return FHT
     */
    double calculate(const RecPmtProp& pmt) {
        // TODO: optimize
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
            fht_1 = m_t_0_1 + m_d_o_wp_1 * constants::inv_c + mag(pmt.pos - m_p_end_wp_1) * constants::inv_c_w - 10.0; // 10 is an arbitrary constant, leads to better performances ==> maybe because of the Tyvek diffusion
            fht_2 = m_t_0_2 + m_d_o_wp_2 * constants::inv_c + mag(pmt.pos - m_p_end_wp_2) * constants::inv_c_w - 10.0;
        }
        return std::min(fht_1, fht_2);
    };

private:

    double m_d_i_wp_1, m_d_o_wp_1;
    vec3 m_p_end_wp_1;

    double m_d_i_wp_2, m_d_o_wp_2;
    vec3 m_p_end_wp_2;

    vec3 m_p_pmt_proj, m_p_1st_light;

    double m_dist_1st_light;
    double m_dist_orig_to_pmt_perp;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_WPFHT_HPP_