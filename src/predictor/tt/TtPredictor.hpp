#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_TT_TTPREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_TT_TTPREDICTOR_HPP_

#include "predictor/Predictor.hpp"

template<ParamsType _Pt>
class TtPredictor : public Predictor<TtMethodTag>, public TrackSetter<_Pt> {

public:

    ~TtPredictor() override = default;

    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        setTrack(params);
        std::transform(first, last, theo, [&](const vec3& hit) { 
            return this->m_orig + dot(hit - this->m_orig, this->m_dir) * this->m_dir;
        });
    }

    ParamsType getIParamsType() override {
        return _Pt;
    }

};

template<>
class TtPredictor<ParamsType::DoubleAcrylic> : public Predictor<TtMethodTag>, public TrackSetter<ParamsType::DoubleAcrylic> {

public:

    ~TtPredictor() override = default;

    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        setTrack(params);
        std::transform(first, last, theo, [&](const vec3& hit) {
            m_hit_1 = m_orig_1 + dot(hit - m_orig_1, m_dir) * m_dir;
            m_hit_2 = m_orig_2 + dot(hit - m_orig_2, m_dir) * m_dir;
            return ( mag2(m_hit_1 - hit) < mag2(m_hit_2 - hit) ) ? m_hit_1 : m_hit_2;
        });
    }

    ParamsType getIParamsType() override {
        return ParamsType::DoubleAcrylic;
    }

protected:

    vec3 m_hit_1, m_hit_2;

};

template<>
class TtPredictor<ParamsType::TripleAcrylic> : public Predictor<TtMethodTag>, public TrackSetter<ParamsType::TripleAcrylic> {

public:

    ~TtPredictor() override = default;

    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        setTrack(params);
        std::transform(first, last, theo, [&](const vec3& hit) {
            m_hit_1 = m_orig_1 + dot(hit - m_orig_1, m_dir) * m_dir;
            m_hit_2 = m_orig_2 + dot(hit - m_orig_2, m_dir) * m_dir;
            m_hit_3 = m_orig_3 + dot(hit - m_orig_3, m_dir) * m_dir;
            return ( mag2(m_hit_1 - hit) < mag2(m_hit_2 - hit) ) ? 
                ( ( mag2(m_hit_1 - hit) < mag2(m_hit_3 - hit) ) ? m_hit_1 : m_hit_3 ) :
                ( ( mag2(m_hit_2 - hit) < mag2(m_hit_3 - hit) ) ? m_hit_2 : m_hit_3 );
        });
    }

    ParamsType getIParamsType() override {
        return ParamsType::TripleAcrylic;
    }

protected:

    vec3 m_hit_1, m_hit_2, m_hit_3;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_TT_TTPREDICTOR_HPP_