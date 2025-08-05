#include "transformer/WpGeomTimeTransformer.hpp"

#include <algorithm>
#include <limits>
#include <new>
#include <unordered_map>

#include "SniperKernel/SniperLog.h"

WpGeomTimeTransformer::WpGeomTimeTransformer(const std::string& name, double time_window_early, double time_window_late, double thold_sep, double thold_q_ratio, double radius_time, double radius_arclength) : 
    Transformer(name, RecPmtType::PMT_WP), 
    m_time_window_early(time_window_early), 
    m_time_window_late(time_window_late), 
    m_thold_sep(thold_sep), 
    m_thold_q_ratio(thold_q_ratio),
    m_radius_time(radius_time), 
    m_radius_arclength(radius_arclength),
    //
    m_dist2_neigh(5000.0 * 5000.0),
    m_neigh_thold(10u),
    m_pmt_thold(3.0),
    m_h_fht(std::make_unique<TH1D>("fht_hist_WpGeom", "FHT Distribution", 200, 0.0, 1000.0)),
    m_width_time_window(35.0),
    //

    // //
    m_qthold_iftime(30.0)
    // //
{
    m_h_fht->SetDirectory(0);
}

void WpGeomTimeTransformer::transform(RecPmtTable& table) {
    if (!findRange(table)) return;

    // ~~~
    // m_nb_neigh.clear();
    // m_nb_neigh.resize(std::distance(m_ftable, m_ltable), 0u);
    // for (RecPmtTable::iterator it = m_ftable; it != m_ltable; ++it) {
    //     if (!it->used || !checkPmtType(*it)) continue;
    //     for (RecPmtTable::iterator jt = m_ftable; jt != m_ltable; ++jt) {
    //         if (!jt->used || !checkPmtType(*jt)) continue;
    //         if (it->id == jt->id) continue;
    //         if (m_dist2_neigh < mag2(it->pos - jt->pos)) continue;
    //         ++m_nb_neigh[std::distance(m_ftable, it)];
    //     }
    // }
    // for (std::size_t k = 0; k < m_nb_neigh.size(); ++k) {
    //     if (m_nb_neigh[k] < m_neigh_thold) m_ftable[k].used = false;
    // }
    // std::cout << "After applying neighbor threshold: " << std::count_if(m_ftable, m_ltable, [&](RecPmtProp& pmt) { return pmt.used; }) << '\n';
    // ~~~

    if (!getEarlyLateTime()) return;
    std::cout << "time_early: " << m_time_early << " time_late: " << m_time_late << '\n';

    for (RecPmtTable::iterator it = m_ftable; it != m_ltable; ++it) {
        if (!it->used || !checkPmtType(*it)) continue;
        if (it->fht < m_time_early || m_time_late + 100.0 < it->fht) it->used = false;
    }

    // ~~~
    // m_first_time_late = m_time_late;
    // for (RecPmtTable::iterator it = m_ftable; it != m_ltable; ++it) {
    //     if (!it->used || !checkPmtType(*it)) continue;
    //     if (it->fht < m_time_early || m_time_late < it->fht) it->used = false;
    // }
    // ~~~

    if (!getEarlyLatePosition()) return;
    std::cout << "pos_early: " << m_pos_early.x << ' ' << m_pos_early.y << ' ' << m_pos_early.z << '\n';
    std::cout << "pos_late: " << m_pos_late.x << ' ' << m_pos_late.y << ' ' << m_pos_late.z << '\n';
    if (!getMaxRange()) return;
    std::cout << "max_range: " << m_max_range << '\n';
    if (!fillClusteredTable()) return;
    // m_thold_nb_neighbors = static_cast<unsigned int>( std::min(0.05 * m_clustered_table.size() + 5.0, 20.0) );
    m_thold_nb_neighbors = static_cast<unsigned int>( std::min(1.0/30.0 * m_clustered_table.size() + 5.0, 15.0) );
    std::cout << "Clustered table size: " << m_clustered_table.size() << " thold_nb_neighbors: " << m_thold_nb_neighbors << '\n';
    dbscan(m_thold_nb_neighbors);
    checkClusters();
    int cluster_id = getClusterWithMaxPmts();
    if (!getMaxChargeInCluster(cluster_id)) return;
    std::cout << "max_q: " << m_max_q << '\n';
    getMeanPosAndTimeInCluster(cluster_id);
    std::cout << "pos_late: " << m_pos_late.x << ' ' << m_pos_late.y << ' ' << m_pos_late.z << '\n';
    std::cout << "time_late: " << m_time_late << '\n';

    m_dist_range = 20050.0 * angle(m_pos_early, m_pos_late);
    m_time_range = m_time_late - m_time_early;

    std::cout << "dist_range: " << m_dist_range << " time_range: " << m_time_range << '\n';

    // Transformer::transform(table);
    std::for_each(table.begin(), table.end(), [&](RecPmtProp& pmt) { transformPmt(pmt); });

    /* std::unordered_map<int, unsigned int> pmt_neigh_cnt;
    for (RecPmtTable::const_iterator it = m_ftable; it != m_ltable; ++it) {
        if (!it->used || !checkPmtType(*it)) continue;
        for (RecPmtTable::const_iterator jt = m_ftable; jt != m_ltable; ++jt) {
            if (!jt->used || !checkPmtType(*jt)) continue;
            if (it->id == jt->id) continue;
            if (mag2(it->pos - it->pos) < m_dist2_neigh) {
                ++pmt_neigh_cnt[it->id];
            }
        }
    }

    for (RecPmtTable::iterator it = m_ftable; it != m_ltable; ++it) {
        if (!it->used || !checkPmtType(*it)) continue;
        if (pmt_neigh_cnt.find(it->id) == pmt_neigh_cnt.end()) it->used = false;
        else if (pmt_neigh_cnt[it->id] <= 2u) it->used = false;
    } */

    table.erase(std::remove_if(m_ftable, m_ltable, [&](RecPmtProp& pmt) { return !pmt.used; }), m_ltable);

    // ~~~
    // std::for_each(m_ftable, m_ltable, [&](RecPmtProp& pmt) { transform(pmt); });
// 
    // m_h_fht->Reset();
    // for (const RecPmtProp& pmt : table) {
    //     if (!pmt.used || !checkPmtType(pmt)) continue;
    //     m_h_fht->Fill(pmt.fht);
    // }
// 
    // int bin_begin = m_h_fht->FindBin(m_first_time_late - m_time_window_late);
    // int bin_end = m_h_fht->FindBin(m_first_time_late);
    // std::cout << "bin_begin: " << bin_begin << " bin_end: " << bin_end << '\n';
    // double time_max_bin = m_first_time_late - m_time_window_late * 0.5;
    // int max_bin = m_h_fht->FindBin(time_max_bin);
    // for (int k = bin_begin; k <= bin_end; ++k) {
    //     if (m_h_fht->GetBinContent(max_bin) < m_h_fht->GetBinContent(k)) {
    //         max_bin = k;
    //     }
    // }
    // time_max_bin = m_h_fht->GetBinCenter(max_bin);
    // std::cout << max_bin << ' ' << time_max_bin << '\n';
    // for (RecPmtProp& pmt : table) {
    //     if (!pmt.used || !checkPmtType(pmt)) continue;
    //     double dist_pmt = arcLengthRatio(pmt, m_dist_range);
    //     if (0.7 <= dist_pmt && (pmt.fht <= time_max_bin - m_width_time_window || time_max_bin + m_width_time_window <= pmt.fht)) pmt.used = false;
    // }
// 
    // std::size_t isize = table.size();
    // table.erase(std::remove_if(m_ftable, m_ltable, [&](RecPmtProp& pmt) { return !pmt.used; }), m_ltable);
    // std::size_t fsize = table.size();
// 
    // LogDebug << isize << " -> " << fsize << " = " << isize - fsize << " PMTs are removed\n";
    // ~~~
}

bool WpGeomTimeTransformer::getEarlyLateTime() {
    m_time_early = std::numeric_limits<double>::infinity();
    m_time_late = -std::numeric_limits<double>::infinity();

    bool found_early = false;
    bool found_late = false;

    for (RecPmtTable::const_iterator it = m_ftable; it != m_ltable; ++it) {
        if (!it->used || !checkPmtType(*it)) continue;
        if (it->q < m_qthold_iftime) continue;
        m_time_early = std::min(m_time_early, it->fht);
        m_time_late = std::max(m_time_late, it->fht);
        found_early = true;
        found_late = true;
    }

    if (!found_early || !found_late) {
        LogWarn << "Early/Late time not found, cannot continue WpGeomTimeTransformer\n";
        return false;
    }

    return true;
}

bool WpGeomTimeTransformer::getEarlyLatePosition() {
    m_pos_early = vec3{0.0, 0.0, 0.0};
    m_pos_late = vec3{0.0, 0.0, 0.0};
    double totq_early = 0.0, totq_late = 0.0;

    bool found_early = false;
    bool found_late = false;

    for (RecPmtTable::const_iterator it = m_ftable; it != m_ltable; ++it) {
        if (!it->used || !checkPmtType(*it)) continue;
        if (m_time_early <= it->fht && it->fht <= m_time_early + m_time_window_early) {
            m_pos_early += it->pos * it->q;
            totq_early += it->q;
            found_early = true;
        }
        if (m_time_late - m_time_window_late <= it->fht && it->fht <= m_time_late) {
            m_pos_late += it->pos * it->q;
            totq_late += it->q;
            found_late = true;
        }
    }

    if (!found_early || !found_late) {
        LogWarn << "Early/Late position not found, cannot continue WpGeomTimeTransformer\n";
        return false;
    }

    m_pos_early /= totq_early;
    m_pos_late /= totq_late;

    if (mag(m_pos_late - m_pos_early) < 5000.0) {
        LogInfo << "Early/Late position are too close, event might be a stopping muon, cannot continue WpGeomTimeTransformer\n";
        return false;
    }

    return true;
}

bool WpGeomTimeTransformer::getMaxRange() {
    m_plan_z = unit(cross(m_pos_late, m_pos_early));

    double min_dist = std::numeric_limits<double>::infinity();
    double max_dist = -std::numeric_limits<double>::infinity();

    bool found_range = false;

    for (RecPmtTable::const_iterator it = m_ftable; it != m_ltable; ++it) {
        if (!it->used || !checkPmtType(*it)) continue;
        double dist_pmt = arcLengthRatio(*it, 1.0);
        min_dist = std::min(min_dist, dist_pmt);
        max_dist = std::max(max_dist, dist_pmt);
        found_range = true;
    }

    if (!found_range) {
        LogWarn << "Range not found, cannot continue WpGeomTimeTransformer\n";
        return false;
    }

    m_max_range = max_dist - min_dist;
    
    return true;
}

bool WpGeomTimeTransformer::fillClusteredTable() {
    m_clustered_table.clear();
    for (RecPmtTable::const_iterator it = m_ftable; it != m_ltable; ++it) {
        if (!it->used || !checkPmtType(*it)) continue;
        double dist_pmt = arcLengthRatio(*it, m_max_range);
        if (m_thold_sep <= dist_pmt) {
            m_clustered_table.push_back({it, 0});
        }
    }
    if (m_clustered_table.empty()) {
        LogWarn << "No pmts to cluster, cannot continue WpGeomTimeTransformer\n";
        return false;
    }
    return true;
}

void WpGeomTimeTransformer::dbscan(unsigned int thold_nb_neighbors) {
    int cur_cluster_id = 0;
    for (std::size_t k = 0; k < m_clustered_table.size(); ++k) {
        if (m_clustered_table[k].cluster_id) continue; // != 0 ==> already visited
        m_clustered_table[k].cluster_id = -1;
        std::vector<std::size_t> neighbors = regionQuery(k);
        if (neighbors.size() < thold_nb_neighbors) continue;
        m_clustered_table[k].cluster_id = ++cur_cluster_id;
        expandCluster(neighbors, cur_cluster_id, thold_nb_neighbors);
    }
}

void WpGeomTimeTransformer::expandCluster(std::vector<std::size_t>& neighbors, int cluster_id, unsigned int thold_nb_neighbors) {
    while (!neighbors.empty()) {
        std::size_t k = neighbors.back();
        neighbors.pop_back();
        if (m_clustered_table[k].cluster_id) continue; // != 0 ==> already visited
        m_clustered_table[k].cluster_id = cluster_id;
        std::vector<std::size_t> new_neighbors = regionQuery(k);
        if (new_neighbors.size() < thold_nb_neighbors) continue;
        for (std::size_t i = 0; i < new_neighbors.size(); ++i) {
            if (std::find(neighbors.begin(), neighbors.end(), new_neighbors[i]) == neighbors.end()) {
                neighbors.push_back(new_neighbors[i]);
            }
        }
    }
}

std::vector<std::size_t> WpGeomTimeTransformer::regionQuery(const std::size_t idx) {
    std::vector<std::size_t> neighbors;
    double cur_dist_pmt = arcLengthRatio(*m_clustered_table[idx].it_pmt, m_max_range);
    for (std::size_t k = 0; k < m_clustered_table.size(); ++k) {
        if (k == idx) continue;
        double dist_pmt = arcLengthRatio(*m_clustered_table[k].it_pmt, m_max_range);
        if (
            (m_clustered_table[idx].it_pmt->fht - m_radius_time <= m_clustered_table[k].it_pmt->fht && m_clustered_table[k].it_pmt->fht <= m_clustered_table[idx].it_pmt->fht + m_radius_time) &&
            (cur_dist_pmt - m_radius_arclength <= dist_pmt && dist_pmt <= cur_dist_pmt + m_radius_arclength)
        ) {
            neighbors.push_back(k);
        }
    }
    return neighbors;
}

void WpGeomTimeTransformer::checkClusters() {
    bool found_cluster = false;
    for (const ClusteredPmt& pmt : m_clustered_table) {
        if (pmt.cluster_id && pmt.cluster_id != -1) {
            found_cluster = true;
            break;
        }
    }
    if (!found_cluster) {
        LogWarn << "No cluster found, trying with thold_nb_neighbor - 10\n";
        for (ClusteredPmt& pmt : m_clustered_table) {
            pmt.cluster_id = 0;
        }
        dbscan( std::max(0u, m_thold_nb_neighbors - 10) );
    }
}

int WpGeomTimeTransformer::getClusterWithMaxPmts() {
    int cluster_with_max_pmts = 0;
    std::unordered_map<int, unsigned int> cluster_counter;
    for (const ClusteredPmt& pmt : m_clustered_table) {
        if (pmt.cluster_id == -1) continue;
        if (cluster_counter.find(pmt.cluster_id) == cluster_counter.end()) {
            cluster_counter[pmt.cluster_id] = 1;
        } 
        else {
            cluster_counter[pmt.cluster_id]++;
        }
    }

    for (const auto& [cluster_id, count] : cluster_counter) {
        if (cluster_counter[cluster_with_max_pmts] < count) {
            cluster_with_max_pmts = cluster_id;
        }
    }
    return cluster_with_max_pmts;
}

bool WpGeomTimeTransformer::getMaxChargeInCluster(int cluster_id) {
    m_max_q = 0.0;
    bool found_maxq = false;
    for (const ClusteredPmt& pmt : m_clustered_table) {
        if (pmt.cluster_id != cluster_id) continue;
        if (m_max_q < pmt.it_pmt->q) {
            m_max_q = pmt.it_pmt->q;
            found_maxq = true;
        }
    }
    if (!found_maxq) {
        LogWarn << "Max charge not found, cannot continue WpGeomTimeTransformer\n";
        return false;
    }
    return true;
}

void WpGeomTimeTransformer::getMeanPosAndTimeInCluster(int cluster_id) {
    m_pos_late = vec3{0.0, 0.0, 0.0};
    double totq_late = 0.0;
    m_time_late = 0.0;

    for (const ClusteredPmt& pmt : m_clustered_table) {
        if (pmt.cluster_id != cluster_id) continue;
        if (m_thold_q_ratio * std::log10(m_max_q) < std::log10(pmt.it_pmt->q)) {
            m_pos_late += pmt.it_pmt->pos * pmt.it_pmt->q;
            totq_late += pmt.it_pmt->q;
            m_time_late += pmt.it_pmt->fht * pmt.it_pmt->q;
        }
    }

    m_pos_late /= totq_late;
    m_time_late /= totq_late;
}

double WpGeomTimeTransformer::arcLengthRatio(const RecPmtProp& pmt, double range) {
    vec3 pmt_proj = pmt.pos - dot(pmt.pos, m_plan_z) * m_plan_z;
    return 20050.0 * angle(m_pos_early, pmt_proj) / range;
}

void WpGeomTimeTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    double dist_pmt = arcLengthRatio(pmt, m_dist_range);
    double time_pmt = (pmt.fht - m_time_early) / m_time_range;
    if (
        // (dist_pmt <= 0.7 && m_time_early + 50.0 <= pmt.fht) ||
        // (0.3 <= dist_pmt && time_pmt <= 0.9) || 
        (dist_pmt <= 0.7 && 0.4 <= time_pmt) ||
        (0.3 <= dist_pmt && time_pmt <= 0.7) ||
        (1.25 <= time_pmt) // \in [1.2, 1.3]
        // (0.9 <= time_pmt && time_pmt <= -0.75 * dist_pmt + 1.65) ||
        // (0.5 <= dist_pmt && 10000.0 <= mag(m_pos_late - pmt.pos))
    ) {
        pmt.used = false;
    }
    return;
}