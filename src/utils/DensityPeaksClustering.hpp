#ifndef CDWPTTCHI2RECTOOL_UTILS_DENSITYPEAKSCLUSTERING_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_DENSITYPEAKSCLUSTERING_HPP_

#include "utils/RecPmtProp.hpp"

class DensityPeaksClustering {

public:

    DensityPeaksClustering();
    ~DensityPeaksClustering() = default;

    void setDistance(double dist) { m_dist_2 = dist * dist; };

    void cluster(const RecPmtTable& table);

    const std::vector<double>& getDensities() const { return m_densities; };
    const std::vector<double>& getDistances() const { return m_distances; };

private:

    double m_dist_2;

    std::vector<double> m_densities;
    std::vector<double> m_distances;

    void fillDensities(const RecPmtTable& table);
    void fillDistances(const RecPmtTable& table);

};

#endif // CDWPTTCHI2RECTOOL_UTILS_DENSITYPEAKSCLUSTERING_HPP_