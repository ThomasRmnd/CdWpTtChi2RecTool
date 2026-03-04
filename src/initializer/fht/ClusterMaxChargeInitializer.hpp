#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERMAXCHARGEINITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERMAXCHARGEINITIALIZER_HPP_

#include <memory>

#include <TH1D.h>

#include "initializer/Initializer.hpp"

/**
 * @class ClusterMaxChargeInitializer
 * 
 * @brief Derived class to calculate track parameters initial guess for FHT method
 * This initial guess is designed for LS phase single track 
 */
class ClusterMaxChargeInitializer : public Initializer<FhtMethodTag> {

public:

    enum class DetectorMode {
        CD_ONLY       = 0b01,
        WP_ONLY       = 0b10,
        CDWP_COMBINED = 0b11
    };

    enum class ProjectionMode {
        ACRYLIC_SPHERE,
        CD_SPHERE
    };

    ClusterMaxChargeInitializer(const std::string& name, int nbins_cd, int nbins_wp, double xmin, double xmax, double window_min, double window_max, double qthold, double offset_final_cd, double exclusion_time_wp, double leading_entries_ratio, DetectorMode detmode, ProjectionMode projmode);

    ~ClusterMaxChargeInitializer() override = default;

    void configure(const SniperJSON& config);

    ParamsType getOParamsType() override;

    /**
     * @brief Calculate the track parameters initial guess
     * 
     * @param data experimental data vector
     * 
     * @return Boolean whether the calculation is successfull or not
     */
    bool initiate(const RecPmtTable& table) override;

private:

    std::unique_ptr<TH1D> m_hist_cd;
    std::unique_ptr<TH1D> m_hist_wp;
    double m_window_min;
    double m_window_max;
    double m_qthold;
    double m_offset_final_cd;
    double m_exclusion_time_wp;
    double m_leading_entries_ratio;
    DetectorMode m_detmode;
    ProjectionMode m_projmode;
    double m_proj_radius;

    bool getCdInit(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable, double qthold, double& itime, vec3& ipos, vec3& fpos);
    bool getWpInit(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable, double qthold, vec3& ipos, vec3& fpos);

};

ClusterMaxChargeInitializer::DetectorMode operator&(const ClusterMaxChargeInitializer::DetectorMode& a, const ClusterMaxChargeInitializer::DetectorMode& b);

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERMAXCHARGEINITIALIZER_HPP_