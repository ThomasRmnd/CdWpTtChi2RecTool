#include "estimator/tt/converter/FuzeNeighborConverter.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

#include "Identifier/TtID.h"

DECLARE_TOOL(FuzeNeighborConverter);

bool FuzeNeighborConverter::convert(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) {
    if (!getWallHits(ftable, ltable)) return false;
    if (!cleanWallHits()) return false;
    return convertToHits();
}

bool FuzeNeighborConverter::getWallHits(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) {
    m_wall_hits.clear();
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        Identifier id = TtID::id(it->id);
        int wall_id = TtID::wall_id(id);
        if (it->pos.x < 1.e9) {
            m_wall_hits[wall_id].x.push_back(it->pos.x);
            m_wall_hits[wall_id].z_x = it->pos.z;
        }
        else if (it->pos.y < 1.e9) {
            m_wall_hits[wall_id].y.push_back(it->pos.y);
            m_wall_hits[wall_id].z_y = it->pos.z;
        }
        else {
            LogError << "Wrong PMT position (id: " << id << ")\n";
            return false;
        }
    }
    return true;
}

bool FuzeNeighborConverter::cleanWallHits() {
    m_cleaned_wall_hits.clear();
    for (std::unordered_map<int, WallHits>::iterator it = m_wall_hits.begin(); it != m_wall_hits.end(); ++it) {
        if (it->second.x.empty() || it->second.y.empty()) continue;
        WallHits wh;
        wh.x = cleanAxisHits(it->second.x);
        wh.y = cleanAxisHits(it->second.y);
        wh.z_x = it->second.z_x;
        wh.z_y = it->second.z_y;
        m_cleaned_wall_hits.push_back(wh);
    }
    return true;
}

std::vector<double> FuzeNeighborConverter::cleanAxisHits(std::vector<double>& axis_hits) {
    std::vector<double> clean_axis_hits;
    std::sort(axis_hits.begin(), axis_hits.end());
    unsigned int nb = 0;
    double coord = 0.0;
    for (std::size_t i = 0; i < axis_hits.size(); ++i) {
        coord += axis_hits[i];
        ++nb;
        if ( (i == axis_hits.size() - 1ul) || (axis_hits[i + 1] - axis_hits[i] > 30.0) ) {
            clean_axis_hits.push_back(coord / nb);
            coord = 0.0;
            nb = 0;
        }
    }
    return clean_axis_hits;
}

bool FuzeNeighborConverter::convertToHits() {
    m_hits.clear();
    for (const WallHits& wh : m_cleaned_wall_hits) {
        double z = (wh.z_x + wh.z_y) / 2.;
        for (std::size_t i = 0; i < wh.x.size(); ++i) {
            for (std::size_t j = 0; j < wh.y.size(); ++j) {
                m_hits.push_back(vec3{wh.x[i], wh.y[j], z});
            }
        }
    }
    if (m_hits.empty()) {
        LogWarn << "There are no hits\n";
    }
    return true;
}