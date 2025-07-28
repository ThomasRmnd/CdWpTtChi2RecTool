#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERBUNDLEINITIALIZER_H_
#define CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERBUNDLEINITIALIZER_H_

#include "initializer/Initializer.hpp"
#include "utils/DensityPeaksClustering.hpp"

/**
 * @class ClusterBundleInitializer
 * 
 * @brief Derived class to calculate track parameters initial guess for FHT method
 * This initial guess is designed for LS phase double track 
 */
class ClusterBundleInitializer : public Initializer<FhtMethodTag> {

public:

    ClusterBundleInitializer(const std::string& name, double dens_thold, double dist_thold);

    ~ClusterBundleInitializer() override = default;

    ParamsType getOParamsType() override;

    void configure(const SniperJSON& config);

    /**
     * @brief Calculate the track parameters initial guess
     * 
     * @param data experimental data vector
     * 
     * @return Boolean whether the calculation is successfull or not
     */
    bool initiate(const RecPmtTable& table) override;

private:

    struct Cluster {
        Cluster() = default;
        Cluster(const vec3& pos_, double charge_, double time_) : pos{pos_}, charge{charge_}, time{time_} {};

        vec3 pos;
        double charge;
        double time;
    };

    DensityPeaksClustering m_dpc;

    double m_dens_thold, m_dist_thold;
    
    std::vector<Cluster> m_clusters;

    bool strategy1Cluster(const RecPmtTable& table);
    bool strategy2Cluster(const RecPmtTable& table);
    bool strategy3Cluster(const RecPmtTable& table);
    bool strategy4Cluster(const RecPmtTable& table);
    bool strategy4ClusterSpe(const RecPmtTable& table);

};

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERBUNDLEINITIALIZER_H_