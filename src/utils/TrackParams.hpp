#ifndef CDWPTTCHI2RECTOOL_UTILS_TRACKPARAMS_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_TRACKPARAMS_HPP_

#include <unordered_map>

#include "utils/constants.hpp"
#include "utils/vec3.hpp"

struct ParamsTag {};
struct SingleAcrylicParamsTag : public ParamsTag {};
struct SingleStoppingAcrylicParamsTag : public ParamsTag {};
struct DoubleAcrylicParamsTag : public ParamsTag {};
struct BundleAcrylicParamsTag : public ParamsTag {};
struct SingleCdParamsTag : public ParamsTag {};
struct SingleTtParamsTag : public ParamsTag {};

enum class ParamsType {
    Unknown,
    SingleAcrylic,         // single through-going muon, with entry point contained in the Acrylic sphere
    SingleStoppingAcrylic, // single stopping muon, with entry point contained in the Acrylic sphere
    DoubleAcrylic,         // double muon (through-going or stopping for now), with entry point at the Acrylic sphere
    BundleAcrylic,         // bundle muon (through-going or stopping for now), with entry point at the Acrylic sphere
    SingleCd,              // single through-going muon, with entry point contained in the CD sphere
    SingleTt,              // single muon (through-going or stopping), with entry point at the TT level
};

class IParamsHandler {

public:

    virtual ~IParamsHandler() = default;

    virtual ParamsType getIParamsType() = 0;

};

class OParamsHandler {

public:

    virtual ~OParamsHandler() = default;

    virtual ParamsType getOParamsType() = 0;

};

template<typename _ParamsTag>
struct TrackSetterHelper;

template<typename _ParamsTag>
struct ParamsTraits;

// SingleAcrylic: t_0, theta_i, phi_i, theta_d, phi_d

template<>
struct TrackSetterHelper<SingleAcrylicParamsTag> {

    static void set(const double* params, double& t_0, vec3& orig, vec3& dir, double& length) {
        t_0 = params[0];
        orig = from_spherical(constants::r_acrylic, params[1], params[2]);
        dir = from_spherical(1.0, params[3], params[4]);
        length = -2.0 * dot(orig, dir);
    }

};

template<>
struct ParamsTraits<SingleAcrylicParamsTag> {

    static constexpr ParamsType type = ParamsType::SingleAcrylic;
    static constexpr std::size_t size = 5ul;
    static constexpr double defaults[5] = {300.0, 0.0, 0.0, -constants::pi, 0.0};
    static constexpr double steps[5] = {5.0, 0.1, 0.1, 0.1, 0.1};
    static constexpr std::string_view names[5] = {"t_0", "theta_i", "phi_i", "theta_d", "phi_d"};

};

// SingleStoppingAcrylic: t_0, theta_i, phi_i, theta_d, phi_d, length

template<>
struct TrackSetterHelper<SingleStoppingAcrylicParamsTag> {

    static void set(const double* params, double& t_0, vec3& orig, vec3& dir, double& length) {
        t_0 = params[0];
        orig = from_spherical(constants::r_acrylic, params[1], params[2]);
        dir = from_spherical(1.0, params[3], params[4]);
        length = -2.0 * dot(orig, dir) * params[5];
    }

};

template<>
struct ParamsTraits<SingleStoppingAcrylicParamsTag> {

    static constexpr ParamsType type = ParamsType::SingleStoppingAcrylic;
    static constexpr std::size_t size = 6ul;
    static constexpr double defaults[6] = {300.0, 0.0, 0.0, -constants::pi, 0.0, 1.0};
    static constexpr double steps[6] = {5.0, 0.1, 0.1, 0.1, 0.1, 0.1};
    static constexpr std::string_view names[6] = {"t_0", "theta_i", "phi_i", "theta_d", "phi_d", "length"};

};

// DoubleAcrylic: t_0_1, theta_i_1, phi_i_1, t_0_2, theta_i_2, phi_i_2, theta_d, phi_d

template<>
struct TrackSetterHelper<DoubleAcrylicParamsTag> {

    // first call to setTrack is for the first track, second call is for the second track
    static void set(const double* params, double& t_0, vec3& orig, vec3& dir, double& length) {
        t_0 = params[idx * 3];
        orig = from_spherical(constants::r_acrylic, params[idx * 3 + 1], params[idx * 3 + 2]);
        dir = from_spherical(1.0, params[6], params[7]);
        length = -2.0 * dot(orig, dir);
        ++idx %= 2;
    }

    inline static std::size_t idx = 0ul;

};

template<>
struct ParamsTraits<DoubleAcrylicParamsTag> {

    static constexpr ParamsType type = ParamsType::DoubleAcrylic;
    static constexpr std::size_t size = 8ul;
    static constexpr double defaults[8] = {300.0, 0.0, 0.0, 300.0, 0.0, 0.0, -constants::pi, 0.0};
    static constexpr double steps[8] = {5.0, 0.1, 0.1, 5., 0.1, 0.1, 0.1, 0.1};
    static constexpr std::string_view names[8] = {"t_0_1", "theta_i_1", "phi_i_1", "t_0_2", "theta_i_2", "phi_i_2", "theta_d", "phi_d"};

};

// BundleAcrylic: t_0_1, theta_i_1, phi_i_1, t_0_2, theta_i_2, phi_i_2, t_0_3, theta_i_3, phi_i_3, theta_d, phi_d
// TODO: change it to resizable number of tracks

template<>
struct TrackSetterHelper<BundleAcrylicParamsTag> {

    // first call to setTrack is for the first track, second call is for the second track, third call is for the third track
    static void set(const double* params, double& t_0, vec3& orig, vec3& dir, double& length) {
        t_0 = params[idx * 3];
        orig = from_spherical(constants::r_acrylic, params[idx * 3 + 1], params[idx * 3 + 2]);
        dir = from_spherical(1.0, params[9], params[10]);
        length = -2.0 * dot(orig, dir);
        ++idx %= 3;
    }

    inline static std::size_t idx = 0ul;

};

template<>
struct ParamsTraits<BundleAcrylicParamsTag> {

    static constexpr ParamsType type = ParamsType::BundleAcrylic;
    static constexpr std::size_t size = 8ul;
    static constexpr double defaults[8] = {300.0, 0.0, 0.0, 300.0, 0.0, 0.0, -constants::pi, 0.0};
    static constexpr double steps[8] = {5.0, 0.1, 0.1, 5., 0.1, 0.1, 0.1, 0.1};
    static constexpr std::string_view names[8] = {"t_0_1", "theta_i_1", "phi_i_1", "t_0_2", "theta_i_2", "phi_i_2", "theta_d", "phi_d"};

};

// SingleCd: t_0, theta_i, phi_i, theta_d, phi_d

template<>
struct TrackSetterHelper<SingleCdParamsTag> {

    static void set(const double* params, double& t_0, vec3& orig, vec3& dir, double& length) {
        t_0 = params[0];
        orig = from_spherical(constants::r_cd, params[1], params[2]);
        dir = from_spherical(1.0, params[3], params[4]);
        length = -2.0 * dot(orig, dir);
    }

};

template<>
struct ParamsTraits<SingleCdParamsTag> {

    static constexpr ParamsType type = ParamsType::SingleCd;
    static constexpr std::size_t size = 5ul;
    static constexpr double defaults[5] = {300.0, 0.0, 0.0, -constants::pi, 0.0};
    static constexpr double steps[5] = {5.0, 0.1, 0.1, 0.1, 0.1};
    static constexpr std::string_view names[5] = {"t_0", "theta_i", "phi_i", "theta_d", "phi_d"};

};

// ParamsType::SingleTt: x, y, z, dx, dy, dz
// TODO: change it to: x, y, 22000.0 (?), dx, dy, 1.0

template<>
struct TrackSetterHelper<SingleTtParamsTag> {

    static void set(const double* params, double& t_0, vec3& orig, vec3& dir, double& length) {
        t_0 = 0.0;
        orig = vec3{params[0], params[1], params[2]};
        dir = unit(vec3{params[3], params[4], params[5]});
        length = -2.0 * dot(orig, dir);
    };

};

template<>
struct ParamsTraits<SingleTtParamsTag> {

    static constexpr ParamsType type = ParamsType::SingleTt;
    static constexpr std::size_t size = 6u; // TODO: change it to 4
    static constexpr double defaults[6] = {0.0, 0.0, 0.0, 0.0, 0.0, -1.0};
    static constexpr double steps[6] = {1.0, 1.0, 0.0, 1.0, 1.0, 0.0}; // 0.0 because fixed variable
    static constexpr const std::string_view names[6] = {"x", "y", "z", "dx", "dy", "dz"};

};

template<typename _ParamsTag>
class TrackSetter {

    static_assert(std::is_base_of<ParamsTag, _ParamsTag>::value, "Tag must derive from ParamsTag");

public:

    virtual ~TrackSetter() = default;

    virtual void setTrack(const double* params) {
        TrackSetterHelper<_ParamsTag>::set(params, m_t_0, m_orig, m_dir, m_length);

    }

protected:

    double m_t_0;
    vec3 m_orig, m_dir;
    double m_length;

};

template<>
class TrackSetter<DoubleAcrylicParamsTag> {

public:

    virtual ~TrackSetter() = default;

    virtual void setTrack(const double* params) {
        TrackSetterHelper<DoubleAcrylicParamsTag>::set(params, m_t_0_1, m_orig_1, m_dir, m_length_1);
        TrackSetterHelper<DoubleAcrylicParamsTag>::set(params, m_t_0_2, m_orig_2, m_dir, m_length_2);
    }

protected:

    double m_t_0_1, m_t_0_2;
    vec3 m_orig_1, m_orig_2, m_dir;
    double m_length_1, m_length_2;

};

template<>
class TrackSetter<BundleAcrylicParamsTag> {

public:

    virtual ~TrackSetter() = default;

    virtual void setTrack(const double* params) {
        TrackSetterHelper<BundleAcrylicParamsTag>::set(params, m_t_0_1, m_orig_1, m_dir, m_length_1);
        TrackSetterHelper<BundleAcrylicParamsTag>::set(params, m_t_0_2, m_orig_2, m_dir, m_length_2);
        TrackSetterHelper<BundleAcrylicParamsTag>::set(params, m_t_0_3, m_orig_3, m_dir, m_length_3);
    }

protected:

    double m_t_0_1, m_t_0_2, m_t_0_3;
    vec3 m_orig_1, m_orig_2, m_orig_3, m_dir;
    double m_length_1, m_length_2, m_length_3;

};

inline static const std::unordered_map<ParamsType, std::size_t> g_track_type_to_size = {
    {ParamsType::SingleAcrylic,         ParamsTraits<SingleAcrylicParamsTag>::size},
    {ParamsType::SingleStoppingAcrylic, ParamsTraits<SingleStoppingAcrylicParamsTag>::size},
    {ParamsType::DoubleAcrylic,         ParamsTraits<DoubleAcrylicParamsTag>::size},
    {ParamsType::BundleAcrylic,         ParamsTraits<BundleAcrylicParamsTag>::size},
    {ParamsType::SingleCd,              ParamsTraits<SingleCdParamsTag>::size},
    {ParamsType::SingleTt,              ParamsTraits<SingleTtParamsTag>::size}
};

#endif // CDWPTTCHI2RECTOOL_UTILS_TRACKPARAMS_HPP_