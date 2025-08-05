#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_CDFHT_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_CDFHT_HPP_

#include "utils/constants.hpp"
#include "utils/RecPmtProp.hpp"
#include "utils/TrackParams.hpp"

/**
 * @class CdFht
 * @tparam _ParamsTag Method tag
 * 
 * @brief Base class to calculate expected FHT for CD PMT.
 */
template<typename _ParamsTag>
class CdFht : public TrackSetter<_ParamsTag> {

    static_assert(std::is_base_of<ParamsTag, _ParamsTag>::value, "Tag must derive from ParamsTag");

public:

    virtual ~CdFht() override = default;

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
 * @class NoRefractionLsCdFht
 * @tparam _ParamsTag Method tag
 * 
 * @brief Derived class to calculate expected FHT for CD PMT.
 * 
 * This calculation is designed for LS phase.
 * It do not take into account refraction between LS and water.
 */
template<typename _ParamsTag>
class NoRefractionLsCdFht final : public CdFht<_ParamsTag> {

public:

    ~NoRefractionLsCdFht() final override = default;

    void setTrack(const double* params) final override {
        CdFht<_ParamsTag>::setTrack(params);
        m_half_length = 0.5 * this->m_length;
        m_t_end = this->m_t_0 + this->m_length / constants::c; // ~ this->m_length * constants::inv_c but even though multiplying by inversion is faster, dividing seems to be more precise
        m_p_end = this->m_orig + this->m_length * this->m_dir;
    }

    /**
     * @brief Calculate the expected FHT
     * 
     * @param pmt PMT information
     * 
     * @return FHT
     */
    double calculate(const RecPmtProp& pmt) final override {
        m_dist_orig_to_pmt_perp = m_half_length + dot(this->m_dir, pmt.pos);
        m_dist_1st_light = m_dist_orig_to_pmt_perp - mag(pmt.pos - (this->m_orig + m_dist_orig_to_pmt_perp * this->m_dir)) * constants::inv_tan_cherenkov_ls;

        if (m_dist_1st_light < 0.0) return this->m_t_0 + constants::inv_c_ls * mag(pmt.pos - this->m_orig);        
        else if (this->m_length < m_dist_1st_light) return m_t_end + constants::inv_c_ls * mag(pmt.pos - m_p_end);
        else return this->m_t_0 + constants::one_minus_n_ls_2_div_c * m_dist_1st_light + constants::n_ls_2_div_c * m_dist_orig_to_pmt_perp;
    }

private:

    double m_half_length;
    double m_t_end;
    vec3 m_p_end;

    double m_dist_orig_to_pmt_perp, m_dist_1st_light;

};

/**
 * @class NoRefractionLsCdFht<DoubleAcrylicParamsTag>
 * 
 * @brief Derived class to calculate expected FHT for CD PMT, specialization with DoubleAcrylicParamsTag.
 * 
 * This calculation is designed for LS phase.
 * It do not take into account refraction between LS and water.
 */
template<>
class NoRefractionLsCdFht<DoubleAcrylicParamsTag> final : public CdFht<DoubleAcrylicParamsTag> {

public:

    ~NoRefractionLsCdFht() final override = default;

    void setTrack(const double* params) final override {
        CdFht<DoubleAcrylicParamsTag>::setTrack(params);
        m_half_length_1 = 0.5 * m_length_1;
        m_t_end_1 = m_t_0_1 + m_length_1 / constants::c; // ~ m_length_1 * constants::inv_c but even though multiplying by inversion is faster, dividing seems to be more precise
        m_p_end_1 = m_orig_1 + m_length_1 * m_dir;
        m_half_length_2 = 0.5 * m_length_2;
        m_t_end_2 = m_t_0_2 + m_length_2 / constants::c; // ~ m_length_2 * constants::inv_c but even though multiplying by inversion is faster, dividing seems to be more precise
        m_p_end_2 = m_orig_2 + m_length_2 * m_dir;
    }

    /**
     * @brief Calculate the expected FHT
     * 
     * @param pmt PMT information
     * 
     * @return FHT
     */
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

private:

    double m_half_length_1, m_half_length_2;
    double m_t_end_1, m_t_end_2;
    vec3 m_p_end_1, m_p_end_2;

    double m_fht_1, m_fht_2;
    double m_dist_orig_to_pmt_perp, m_dist_1st_light;

};

/**
 * @class WaterPhaseCdFht
 * @tparam _ParamsTag Method tag
 * 
 * @brief Derived class to calculate expected FHT for CD PMT.
 * 
 * This calculation is designed for water phase.
 */
template<typename _ParamsTag>
class WaterPhaseCdFht final : public CdFht<_ParamsTag> {

public:

    ~WaterPhaseCdFht() final override = default;

    void setTrack(const double* params) final override {
        CdFht<_ParamsTag>::setTrack(params);
        m_half_length = 0.5 * this->m_length;
    }

    /**
     * @brief Calculate the expected FHT
     * 
     * @param pmt PMT information
     * 
     * @return FHT
     */
    double calculate(const RecPmtProp& pmt) final override {
        m_dist_orig_to_pmt_perp = m_half_length + dot(this->m_dir, pmt.pos);
        m_dist_1st_light = m_dist_orig_to_pmt_perp - mag(pmt.pos - (this->m_orig + m_dist_orig_to_pmt_perp * this->m_dir)) * constants::inv_tan_cherenkov_w;
        return this->m_t_0 + constants::one_minus_n_w_2_div_c * m_dist_1st_light + constants::n_w_2_div_c * m_dist_orig_to_pmt_perp;
    }

private:

    double m_half_length;

    double m_dist_orig_to_pmt_perp, m_dist_1st_light;

};


#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_FHT_CDFHT_HPP_