#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_CDFHT_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_CDFHT_HPP_

#include "utils/TrackParams.hpp"

#include "utils/constants.hpp"
#include "utils/RecPmtProp.hpp"

template<ParamsType _Pt>
class CdFht : public TrackSetter<_Pt> {

public:

    virtual ~CdFht() override = default;

    virtual double calculate(const RecPmtProp& pmt) = 0;

};

template<ParamsType _Pt>
class NoRefractionLsCdFht final : public CdFht<_Pt> {

public:

    ~NoRefractionLsCdFht() final override = default;

    double calculate(const RecPmtProp& pmt) final override {
        m_dist_orig_to_pmt_perp = m_half_length + dot(this->m_dir, pmt.pos);
        m_dist_1st_light = m_dist_orig_to_pmt_perp - mag(pmt.pos - (this->m_orig + m_dist_orig_to_pmt_perp * this->m_dir)) * constants::inv_tan_cherenkov_ls;

        if (m_dist_1st_light < 0.0) return this->m_t_0 + constants::inv_c_ls * mag(pmt.pos - this->m_orig);        
        else if (this->m_length < m_dist_1st_light) return m_t_end + constants::inv_c_ls * mag(pmt.pos - m_p_end);
        else return this->m_t_0 + constants::one_minus_n_ls_2_div_c * m_dist_1st_light + constants::n_ls_2_div_c * m_dist_orig_to_pmt_perp;
    }

    void setTrack(const double* params) final override {
        CdFht<_Pt>::setTrack(params);
        m_half_length = 0.5 * this->m_length;
        m_t_end = this->m_t_0 + this->m_length / constants::c; // ~ this->m_length * constants::inv_c but even though multiplying by inversion is much faster, dividing seems to be more precise
        m_p_end = this->m_orig + this->m_length * this->m_dir;
    }

private:

    double m_half_length;
    double m_t_end;
    vec3 m_p_end;

    double m_dist_orig_to_pmt_perp, m_dist_1st_light;

};

template<>
class NoRefractionLsCdFht<ParamsType::DoubleAcrylic> final : public CdFht<ParamsType::DoubleAcrylic> {

public:

    ~NoRefractionLsCdFht() final override = default;

    double calculate(const RecPmtProp& pmt) final override {
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

        return std::min(m_fht_1, m_fht_2);
    }

    void setTrack(const double* params) final override {
        CdFht<ParamsType::DoubleAcrylic>::setTrack(params);
        m_half_length_1 = 0.5 * m_length_1;
        m_t_end_1 = m_t_0_1 + m_length_1 / constants::c; // ~ m_length_1 * constants::inv_c but even though multiplying by inversion is much faster, dividing seems to be more precise
        m_p_end_1 = m_orig_1 + m_length_1 * m_dir;
        m_half_length_2 = 0.5 * m_length_2;
        m_t_end_2 = m_t_0_2 + m_length_2 / constants::c; // ~ m_length_2 * constants::inv_c but even though multiplying by inversion is much faster, dividing seems to be more precise
        m_p_end_2 = m_orig_2 + m_length_2 * m_dir;
    }

private:

    double m_half_length_1, m_half_length_2;
    double m_t_end_1, m_t_end_2;
    vec3 m_p_end_1, m_p_end_2;

    double m_fht_1, m_fht_2;
    double m_dist_orig_to_pmt_perp, m_dist_1st_light;

};

template<>
class NoRefractionLsCdFht<ParamsType::TripleAcrylic> final : public CdFht<ParamsType::TripleAcrylic> {

public:

    ~NoRefractionLsCdFht() final override = default;

    double calculate(const RecPmtProp& pmt) final override {
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

        m_dist_orig_to_pmt_perp = m_half_length_3 + dot(m_dir, pmt.pos);
        m_dist_1st_light = m_dist_orig_to_pmt_perp - mag(pmt.pos - (m_orig_3 + m_dist_orig_to_pmt_perp * m_dir)) * constants::inv_tan_cherenkov_ls;
        if (m_dist_1st_light < 0.0) m_fht_3 = m_t_0_3 + constants::inv_c_ls * mag(pmt.pos - m_orig_3);        
        else if (m_length_3 < m_dist_1st_light) m_fht_3 = m_t_end_3 + constants::inv_c_ls * mag(pmt.pos - m_p_end_3);
        else m_fht_3 = m_t_0_3 + constants::one_minus_n_ls_2_div_c * m_dist_1st_light + constants::n_ls_2_div_c * m_dist_orig_to_pmt_perp;

        return std::min(m_fht_1, m_fht_2, m_fht_3);
    }

    void setTrack(const double* params) final override {
        CdFht<ParamsType::TripleAcrylic>::setTrack(params);
        m_half_length_1 = 0.5 * m_length_1;
        m_t_end_1 = m_t_0_1 + m_length_1 / constants::c; // ~ m_length_1 * constants::inv_c but even though multiplying by inversion is much faster, dividing seems to be more precise
        m_p_end_1 = m_orig_1 + m_length_1 * m_dir;
        m_half_length_2 = 0.5 * m_length_2;
        m_t_end_2 = m_t_0_2 + m_length_2 / constants::c; // ~ m_length_2 * constants::inv_c but even though multiplying by inversion is much faster, dividing seems to be more precise
        m_p_end_2 = m_orig_2 + m_length_2 * m_dir;
        m_half_length_3 = 0.5 * m_length_3;
        m_t_end_3 = m_t_0_3 + m_length_3 / constants::c; // ~ m_length_3 * constants::inv_c but even though multiplying by inversion is much faster, dividing seems to be more precise
        m_p_end_3 = m_orig_3 + m_length_3 * m_dir;
    }

private:

    double m_half_length_1, m_half_length_2, m_half_length_3;
    double m_t_end_1, m_t_end_2, m_t_end_3;
    vec3 m_p_end_1, m_p_end_2, m_p_end_3;

    double m_fht_1, m_fht_2, m_fht_3;
    double m_dist_orig_to_pmt_perp, m_dist_1st_light;

};

template<ParamsType _Pt>
class WaterPhaseCdFht final : public CdFht<_Pt> {

public:

    ~NoRefractionLsCdFht() final override = default;

    double calculate(const RecPmtProp& pmt) final override {
        m_dist_orig_to_pmt_perp = m_half_length + dot(this->m_dir, pmt.pos);
        m_dist_1st_light = m_dist_orig_to_pmt_perp - mag(pmt.pos - (this->m_orig + m_dist_orig_to_pmt_perp * this->m_dir)) * this->m_inv_tan_theta_water;
        return this->m_t_0 + this->m_one_minus_n_water_squared_div_c * m_dist_1st_light + this->m_n_water_squared_div_c * m_dist_orig_to_pmt_perp;
    }

    void setTrack(const double* params) final override {
        CdFht<_Pt>::setTrack(params);
        m_half_length = 0.5 * this->m_length;
    }

private:

    double m_half_length;

    double m_dist_orig_to_pmt_perp, m_dist_1st_light;

};


#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_CDFHT_HPP_