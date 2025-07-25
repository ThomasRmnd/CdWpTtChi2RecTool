#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_TT_TTPREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_TT_TTPREDICTOR_HPP_

#include "predictor/Predictor.hpp"

/**
 * @class TtPredictor
 * @tparam _ParamsTag 
 * 
 * @brief Derived class to calculate expected data for TT method
 */
template<typename _ParamsTag>
class TtPredictor : public Predictor<TtMethodTag>, public TrackSetter<_ParamsTag> {

    static_assert(std::is_base_of<ParamsTag, _ParamsTag>::value, "Tag must derive from ParamsTag");

public:

    ~TtPredictor() override = default;

    ParamsType getIParamsType() override {
        return ParamsTraits<_ParamsTag>::type;
    }

    /**
     * @brief Calculate the expected data
     * 
     * @param first the beginning iterator of the experimental data 
     * @param last the end iterator of the experimental data
     * @param theo the beginning iterator of the expected data
     * @param params track parameters
     * 
     * @return Raw \f$ \chi^2 \f$
     */
    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        this->setTrack(params);
        std::transform(first, last, theo, [&](const vec3& hit) { 
            return this->m_orig + dot(hit - this->m_orig, this->m_dir) * this->m_dir;
        });
    }

};

/**
 * @class TtPredictor<DoubleAcrylicParamsTag>
 * 
 * @brief Derived class for chost calculation for the TT method, specialization with DoubleAcrylicParamsTag
 */
template<>
class TtPredictor<DoubleAcrylicParamsTag> : public Predictor<TtMethodTag>, public TrackSetter<DoubleAcrylicParamsTag> {

public:

    ~TtPredictor() override = default;

    ParamsType getIParamsType() override {
        return ParamsType::DoubleAcrylic;
    }

    /**
     * @brief Calculate the expected data
     * 
     * @param first the beginning iterator of the experimental data 
     * @param last the end iterator of the experimental data
     * @param theo the beginning iterator of the expected data
     * @param params track parameters
     * 
     * @return Raw \f$ \chi^2 \f$
     */
    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        setTrack(params);
        std::transform(first, last, theo, [&](const vec3& hit) {
            m_hit_1 = m_orig_1 + dot(hit - m_orig_1, m_dir) * m_dir;
            m_hit_2 = m_orig_2 + dot(hit - m_orig_2, m_dir) * m_dir;
            return ( mag2(m_hit_1 - hit) < mag2(m_hit_2 - hit) ) ? m_hit_1 : m_hit_2;
        });
    }

protected:

    vec3 m_hit_1, m_hit_2;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_TT_TTPREDICTOR_HPP_