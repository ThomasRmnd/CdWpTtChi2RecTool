#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRPARAM_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRPARAM_HPP_

#include "utils/RecPmtProp.hpp"
#include "utils/vec3.hpp"

class CorrParam {

public:

    virtual double operator()(const RecPmtProp& pmt) = 0;

    virtual void setTrack(const vec3& orig, const vec3& dir) {
        m_orig = orig;
        m_dir = dir;
    };

protected:

    vec3 m_orig; 
    vec3 m_dir; 

};

// angle formed by the PMT, the center of the CD and the track, in the plane perpendicular to the track
class AngleCorrParam final : public CorrParam {

public:

    double operator()(const RecPmtProp& pmt) final override {
        m_cntr_to_pmt = pmt.pos - dot(pmt.pos, m_dir) * m_dir;
        return angle(m_cntr_to_pmt, m_cntr_to_trk);
    };

    void setTrack(const vec3& orig, const vec3& dir) final override {
        CorrParam::setTrack(orig, dir);
        m_cntr_to_trk = m_orig - dot(m_orig, m_dir) * m_dir;
    };

protected:

    vec3 m_cntr_to_trk;
    vec3 m_cntr_to_pmt;

};

// cosine of the angle defined just above
class CosAngleCorrParam final : public CorrParam {

public:

    double operator()(const RecPmtProp& pmt) final override {
        m_cntr_to_pmt = pmt.pos - dot(pmt.pos, m_dir) * m_dir;
        return dot(m_cntr_to_pmt, m_cntr_to_trk) / mag(m_cntr_to_pmt);
    };

    void setTrack(const vec3& orig, const vec3& dir) final override {
        CorrParam::setTrack(orig, dir);
        m_cntr_to_trk = m_orig - dot(m_orig, m_dir) * m_dir;
        m_cntr_to_trk = unit(m_cntr_to_trk);
    };

protected:

    vec3 m_cntr_to_trk;
    vec3 m_cntr_to_pmt;

};

// distance between the track init point and the projection of the PMT onto the track
class DistProjPmtToOrigCorrParam final : public CorrParam {

public:

    double operator()(const RecPmtProp& pmt) final override {
        return dot(pmt.pos, m_dir) + m_half_length;
    };

    void setTrack(const vec3& orig, const vec3& dir) final override {
        CorrParam::setTrack(orig, dir);
        m_half_length = - dot(m_orig, m_dir);
    };

protected:

    double m_half_length;

};

// minimum distance between the track and the center of the CD
class DistTrackToCenterCorrParam final : public CorrParam {

public:

    double operator()(const RecPmtProp& pmt) final override {
        (void)pmt;
        return m_dist_to_cntr;
    };

    void setTrack(const vec3& orig, const vec3& dir) final override {
        CorrParam::setTrack(orig, dir);
        m_dist_to_cntr = mag(cross(m_orig, m_dir));
    };

protected:

    double m_dist_to_cntr;

};


// squared of the distance defined just above
class DistTrackToCenterSquaredCorrParam final : public CorrParam {

public:

    double operator()(const RecPmtProp& pmt) final override {
        (void)pmt;
        return m_dist_to_cntr_squared;
    };

    void setTrack(const vec3& orig, const vec3& dir) final override {
        CorrParam::setTrack(orig, dir);
        m_dist_to_cntr_squared = mag2(cross(m_orig, m_dir));
    };

protected:

    double m_dist_to_cntr_squared;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRPARAM_HPP_