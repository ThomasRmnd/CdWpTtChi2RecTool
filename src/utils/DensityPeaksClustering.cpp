#include "utils/DensityPeaksClustering.hpp"

DensityPeaksClustering::DensityPeaksClustering() : m_dist_2(3000.0 * 3000.0) {}

void DensityPeaksClustering::cluster(const RecPmtTable& table) {
    fillDensities(table);
    fillDistances(table);
}

void DensityPeaksClustering::fillDensities(const RecPmtTable& table) {
    m_densities.clear();
    m_densities.reserve(table.size());
    for (RecPmtTable::const_iterator it = table.begin(); it != table.end(); ++it) {
        if (!it->used) {
            m_densities.push_back(0.0);
            continue;
        }
        double density = 0.0;
        unsigned int count = 0u;
        for (RecPmtTable::const_iterator jt = table.begin(); jt != table.end(); ++jt) {
            if (!jt->used) continue;
            double dist_2 = mag2(it->pos - jt->pos);
            if (dist_2 < m_dist_2) {
                density += jt->q;
                ++count;
            }
        }
        if (!count) m_densities.push_back(0.0);
        else m_densities.push_back(density / count);
    }
}

void DensityPeaksClustering::fillDistances(const RecPmtTable& table) {
    m_distances.clear();
    m_distances.reserve(table.size());
    for (RecPmtTable::const_iterator it = table.begin(); it != table.end(); ++it) {
        if (!it->used) {
            m_distances.push_back(0.0);
            continue;
        }
        std::size_t i = std::distance(table.begin(), it);
        double min_dist_2 = std::numeric_limits<double>::infinity();
        double max_dist_2 = 0.0;
        bool found = false, found_max = false;
        for (RecPmtTable::const_iterator jt = table.begin(); jt != table.end(); ++jt) {
            if (!jt->used) continue;
            if (it == jt) continue;
            double dist_2 = mag2(it->pos - jt->pos);
            if (max_dist_2 < dist_2) {
                max_dist_2 = dist_2;
                found_max = true;
            }
            std::size_t j = std::distance(table.begin(), jt);
            if (m_densities[j] < m_densities[i]) continue;
            if (dist_2 < min_dist_2) {
                min_dist_2 = dist_2;
                found = true;
            }
        }
        if (!found) {
            if (found_max) m_distances.push_back(std::sqrt(max_dist_2));
            else m_distances.push_back(0.0);
        }
        else m_distances.push_back(std::sqrt(min_dist_2));
    }
}