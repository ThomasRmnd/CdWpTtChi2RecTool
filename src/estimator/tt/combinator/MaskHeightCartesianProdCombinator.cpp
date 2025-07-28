#include "estimator/tt/combinator/MaskHeightCartesianProdCombinator.hpp"

#include "SniperKernel/SniperLog.h"

MaskHeightCartesianProdCombinator::MaskHeightCartesianProdCombinator(const std::string& name, std::size_t max_nb_heights) :
    Combinator(name),
    m_max_nb_heights(max_nb_heights)
{};

bool MaskHeightCartesianProdCombinator::combine(const std::vector<vec3>& hits) {
    m_hits_comb.clear();
    if (!getHeightMap(hits)) return false;
    
    std::vector<vec3> curr;
    std::size_t max_pts = std::min(m_max_nb_heights, m_height_map.size());
    for (std::size_t nb_pts = max_pts; nb_pts > 2; --nb_pts) {
        m_mask.assign(m_height_map.size(), false);
        std::fill(m_mask.end() - nb_pts, m_mask.end(), true);
        do {
            LogDebug << "Mask: ";
            for (bool b : m_mask) {
                std::cout << b << ' ';
            }
            std::cout << '\n';
            if (!cartesianProduct(curr, m_height_map.begin(), m_mask.begin())) return false;
        } while (std::next_permutation(m_mask.begin(), m_mask.end()));
    }
    return true;
}

bool MaskHeightCartesianProdCombinator::getHeightMap(const std::vector<vec3>& hits) {
    m_height_map.clear();
    for (const vec3& h : hits) {
        m_height_map[h.z].push_back(h);
    }
    if (m_height_map.empty()) {
        LogWarn << "The height map is empty\n";
    }
    LogDebug << "Height map: \n";
    for (const auto& it : m_height_map) {
        LogDebug << it.first << " -> ";
        for (const vec3& h : it.second) {
            std::cout << '(' << h.x << ", " << h.y << ", " << h.z << "), ";
        }
        std::cout << '\n';
    }
    return true;
}

bool MaskHeightCartesianProdCombinator::cartesianProduct(std::vector<vec3>& curr, std::unordered_map<double, std::vector<vec3>>::const_iterator it, std::vector<bool>::const_iterator mask_it) {
    if (it == m_height_map.end()) {
        m_hits_comb.push_back(curr);
        if (curr.empty()) {
            LogWarn << "The current combination is empty\n";
        }
        /* LogDebug << "Current combination: ";
        for (const vec3& h : curr) {
            std::cout << '(' << h.x << ", " << h.y << ", " << h.z << "), ";
        }
        std::cout << '\n'; */
        return true;
    }
    std::unordered_map<double, std::vector<vec3>>::const_iterator it_next = it;
    ++it_next;
    if (!(*mask_it)) {
        if (!cartesianProduct(curr, it_next, mask_it + 1)) return false;
        return true;
    }
    for (const vec3& h :it->second) {
        curr.push_back(h);
        cartesianProduct(curr, it_next, mask_it + 1);
        curr.pop_back();
    }
    return true;
}