#include "initializer/fht/ClusterBundleInitializer.hpp"

#include <algorithm>

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(ClusterBundleInitializer);

ClusterBundleInitializer::ClusterBundleInitializer(const std::string& name) : 
    Initializer<FhtMethodTag>(name),
    m_dpc{},
    m_dens_thold{1000.0},
    m_dist_thold{5000.0}
{}

ClusterBundleInitializer::ClusterBundleInitializer(const std::string& name, double dens_thold, double dist_thold) : 
    Initializer<FhtMethodTag>(name),
    m_dpc{},
    m_dens_thold{dens_thold},
    m_dist_thold{dist_thold}
{}

bool ClusterBundleInitializer::initiate(const RecPmtTable& table) {
    unsigned int count = std::count_if(table.begin(), table.end(), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_20INCH); });
    LogDebug << count << " PMTs are used for the initialization\n";

    RecPmtTable dpc_table;
    dpc_table.reserve(count);
    for (const RecPmtProp& pmt : table) {
        if (!hasPmtType(pmt, RecPmtType::PMT_20INCH)) continue;
        dpc_table.push_back(pmt);
    }

    m_dpc.cluster(dpc_table);

    const std::vector<double>& densities = m_dpc.getDensities();
    const std::vector<double>& distances = m_dpc.getDistances();

    m_clusters.clear();

    for (std::size_t k = 0; k < densities.size(); ++k) {
        if (densities[k] < m_dens_thold) continue;
        if (distances[k] < m_dist_thold) continue;
        m_clusters.emplace_back(dpc_table[k].pos, dpc_table[k].totq, dpc_table[k].fht);
    }

    if (m_clusters.size() == 1) return strategy1Cluster(table);
    else if (m_clusters.size() == 2) return strategy2Cluster(table);
    else if (m_clusters.size() == 3) return strategy3Cluster(table);
    else if (m_clusters.size() == 4) return strategy4Cluster(table);
    else {
        LogWarn << "Wrong number of clusters: " << m_clusters.size() << ", skip this event\n";
        return false;
    }
}

bool ClusterBundleInitializer::strategy1Cluster(const RecPmtTable& table) {
    LogDebug << "Strategy 1: 1 cluster\n";

    return false; // TOOD: remove this line

    /* m_dpc.setDistance(1500.0);

    RecPmtTable dpc_table;
    dpc_table.reserve(17612u);
    for (const RecPmtProp& pmt : table) {
        if (!check(pmt)) continue;
        dpc_table.push_back(pmt);
    }

    auto [pmt_qmin, pmt_qmax] = std::minmax_element(dpc_table.begin(), dpc_table.end(), [](const RecPmtProp& lhs, const RecPmtProp& rhs) { return lhs.totq < rhs.totq; });
    double q_thold = std::sqrt(pmt_qmin->q * pmt_qmax->q); // geometric mean
    dpc_table.erase(std::remove_if(dpc_table.begin(), dpc_table.end(), [q_thold](const RecPmtProp& pmt) { return pmt.totq < q_thold; }), dpc_table.end());

    m_dpc.cluster(dpc_table);
    
    const std::vector<double>& densities = m_dpc.getDensities();
    const std::vector<double>& distances = m_dpc.getDistances();

    m_dpc.setDistance(3000.0); // reset distance
    
    std::vector<unsigned int> nb_neighbors(dpc_table.size(), 0u);
    for (std::size_t i = 0; i < dpc_table.size(); ++i) {
        for (std::size_t j = 0; j < dpc_table.size(); ++j) {
            if (i == j) continue;
            if (mag(dpc_table[i].pos - dpc_table[j].pos) < 3000.0) {
                ++nb_neighbors[i];
            }
        }
    }

    m_clusters.clear();
    for (std::size_t k = 0; k < densities.size(); ++k) {
        if (nb_neighbors[k] < 80u) continue;
        if (distances[k] < 1000.0) continue;
        m_clusters.emplace_back(dpc_table[k].pos, densities[k], dpc_table[k].fht);
    }

    std::sort(m_clusters.begin(), m_clusters.end(), [](const Cluster& lhs, const Cluster& rhs) { return lhs.charge > rhs.charge; }); // descending order

    std::vector<bool> notvalid(m_clusters.size(), false);
    for (std::size_t i = 0; i < m_clusters.size(); ++i) {
        if (notvalid[i]) continue;
        for (std::size_t j = i + 1; j < m_clusters.size(); ++j) {
            if (notvalid[j]) continue;
            if (mag(m_clusters[i].pos - m_clusters[j].pos) < 4000.0) {
                notvalid[j] = true;
            }
        }
    }

    LogDebug << "Number of clusters: " << m_clusters.size() << '\n';
    m_clusters.erase(std::remove_if(m_clusters.begin(), m_clusters.end(), [&](const Cluster& cluster) { return notvalid[&cluster - &m_clusters[0]]; }), m_clusters.end());
    LogDebug << "Number of valid clusters: " << m_clusters.size() << '\n';

    if (m_clusters.size() == 2) return strategy2Cluster(table);
    else if (m_clusters.size() == 3) return strategy3Cluster(table);
    else if (m_clusters.size() >= 4) {
        // Keep only the 4 clusters with the highest charge
        m_clusters.resize(4);
        return strategy4ClusterSpe(table);
    }
    else if (m_clusters.size() == 0) {
        LogWarn << "No cluster found, skip this event\n";
        return false;
    }
    // else there is only 1 cluster
    

    if (!m_cmc(table)) return false;
    std::vector<double> params = m_cmc.getIVars();
    vec3 ipos, dir;
    ipos = from_spherical(17700.0, params[1], params[2]);
    dir = from_spherical(1.0, params[3], params[4]);
    vec3 cpos = ipos - dot(ipos, dir) * dir;
    vec3 xaxis = unit(cpos);
    vec3 zaxis = -dir;
    vec3 yaxis = cross(zaxis, xaxis);
    vec3 ipos1 = ipos + 1500.0 * yaxis;
    vec3 ipos2 = ipos - 1500.0 * yaxis;
    m_ivars = std::vector<double>{
        params[0], theta(ipos1), phi(ipos1), 
        params[0], theta(ipos2), phi(ipos2),
        theta(dir), phi(dir)
    };
    return true; */
}

bool ClusterBundleInitializer::strategy2Cluster(const RecPmtTable& table) {
    LogDebug << "Strategy 2: 2 clusters\n";
    // TODO: implement strategy 2
    (void) table;
    return false;
}

bool ClusterBundleInitializer::strategy3Cluster(const RecPmtTable& table) {
    LogDebug << "Strategy 3: 3 clusters\n";
    // TODO: implement strategy 3
    (void) table;

    return false; // TOOD: remove this line

    /* std::sort(m_clusters.begin(), m_clusters.end(), [](const Cluster& lhs, const Cluster& rhs) { return lhs.time < rhs.time; });
    vec3 ipos1 = m_clusters[0].pos;
    vec3 ipos2 = m_clusters[1].pos;
    vec3 fpos = m_clusters[2].pos;

    std::cout << "Cluster 1: (" << ipos1.X() << ", " << ipos1.Y() << ", " << ipos1.Z() << "), " << m_clusters[0].charge << ", " << m_clusters[0].time << '\n';
    std::cout << "Cluster 2: (" << ipos2.X() << ", " << ipos2.Y() << ", " << ipos2.Z() << "), " << m_clusters[1].charge << ", " << m_clusters[1].time << '\n';
    std::cout << "Cluster 3: (" << fpos.X() << ", " << fpos.Y() << ", " << fpos.Z() << "), " << m_clusters[2].charge << ", " << m_clusters[2].time << '\n';

    // TODO: check if the following code is correct
    double dt1 = m_clusters[2].time - m_clusters[0].time;
    double dt2 = m_clusters[2].time - m_clusters[1].time;

    double dist1 = mag(ipos1 - fpos);
    double dist2 = mag(ipos2 - fpos);

    double diff1 = std::abs( dist1 / 299.792458 - dt1 );
    double diff2 = std::abs( dist2 / 299.792458 - dt2 );

    if (diff2 < diff1) std::swap(ipos1, ipos2);

    vec3 dir = unit(fpos - ipos1);

    m_ivars = std::vector<double>{
        m_clusters[0].time, theta(ipos1), phi(ipos1),
        m_clusters[1].time, theta(ipos2), phi(ipos2),
        theta(dir), phi(dir)
    };

    return true; */
}

bool ClusterBundleInitializer::strategy4Cluster(const RecPmtTable& table) {
    LogDebug << "Strategy 4: 4 clusters\n";
    (void) table;

    std::sort(m_clusters.begin(), m_clusters.end(), [](const Cluster& lhs, const Cluster& rhs) { return lhs.time < rhs.time; });
    vec3 ipos1 = m_clusters[0].pos;
    vec3 ipos2 = m_clusters[1].pos;
    vec3 fpos1 = m_clusters[2].pos;
    vec3 fpos2 = m_clusters[3].pos;

    std::cout << "Cluster 1: (" << ipos1.x << ", " << ipos1.y << ", " << ipos1.z << "), " << m_clusters[0].charge << ", " << m_clusters[0].time << '\n';
    std::cout << "Cluster 2: (" << ipos2.x << ", " << ipos2.y << ", " << ipos2.z << "), " << m_clusters[1].charge << ", " << m_clusters[1].time << '\n';
    std::cout << "Cluster 3: (" << fpos1.x << ", " << fpos1.y << ", " << fpos1.z << "), " << m_clusters[2].charge << ", " << m_clusters[2].time << '\n';
    std::cout << "Cluster 4: (" << fpos2.x << ", " << fpos2.y << ", " << fpos2.z << "), " << m_clusters[3].charge << ", " << m_clusters[3].time << '\n';

    double angle1 = angle(fpos1 - ipos1, fpos2 - ipos2);
    double angle2 = angle(fpos2 - ipos1, fpos1 - ipos2);

    std::cout << "Angle 1: " << angle1 * 180.0 / M_PI << ", Angle 2: " << angle2 * 180.0 / M_PI << '\n';

    if (angle2 < angle1) std::swap(fpos1, fpos2);

    vec3 dir = ( (fpos1 - ipos1) + (fpos2 - ipos2) ) / 2.0 ;

    m_params = std::vector<double>{
        m_clusters[0].time, theta(ipos1), phi(ipos1),
        m_clusters[1].time, theta(ipos2), phi(ipos2),
        theta(dir), phi(dir)
    };

    return true;
}

bool ClusterBundleInitializer::strategy4ClusterSpe(const RecPmtTable& table) {
    LogDebug << "Strategy 4: 4 clusters\n";
    (void) table;

    std::sort(m_clusters.begin(), m_clusters.end(), [](const Cluster& lhs, const Cluster& rhs) { return lhs.time < rhs.time; });
    vec3 ipos1 = m_clusters[0].pos;
    vec3 ipos2 = m_clusters[1].pos;
    vec3 fpos1 = m_clusters[2].pos;
    vec3 fpos2 = m_clusters[3].pos;

    std::cout << "Cluster 1: (" << ipos1.x << ", " << ipos1.y << ", " << ipos1.z << "), " << m_clusters[0].charge << ", " << m_clusters[0].time << '\n';
    std::cout << "Cluster 2: (" << ipos2.x << ", " << ipos2.y << ", " << ipos2.z << "), " << m_clusters[1].charge << ", " << m_clusters[1].time << '\n';
    std::cout << "Cluster 3: (" << fpos1.x << ", " << fpos1.y << ", " << fpos1.z << "), " << m_clusters[2].charge << ", " << m_clusters[2].time << '\n';
    std::cout << "Cluster 4: (" << fpos2.x << ", " << fpos2.y << ", " << fpos2.z << "), " << m_clusters[3].charge << ", " << m_clusters[3].time << '\n';

    double angle1 = angle(fpos1 - ipos1, fpos2 - ipos2);
    double angle2 = angle(fpos2 - ipos1, fpos1 - ipos2);

    std::cout << "Angle 1: " << angle1 * 180.0 / M_PI << ", Angle 2: " << angle2 * 180.0 / M_PI << '\n';

    if (angle2 < angle1) std::swap(fpos1, fpos2);

    vec3 dir = ( (fpos1 - ipos1) + (fpos2 - ipos2) ) / 2.0 ;

    m_params = std::vector<double>{
        m_clusters[0].time, theta(ipos1), phi(ipos1),
        m_clusters[1].time, theta(ipos2), phi(ipos2),
        theta(dir), phi(dir)
    };

    return true;
}

ParamsType ClusterBundleInitializer::getOParamsType() {
    return ParamsType::DoubleAcrylic;
}