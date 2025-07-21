#include "transformer/TtCrossTalkTransformer.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

#include "Identifier/TtID.h"

DECLARE_TOOL(TtCrossTalkTransformer);

TtCrossTalkTransformer::TtCrossTalkTransformer(const std::string& name) : 
    Transformer(name, RecPmtType::PMT_TT)
{}

void TtCrossTalkTransformer::transform(RecPmtTable& table) {
    if (!findRange(table)) return;
    getHitsStrips();
    filterCrossTalk();
    
    std::size_t isize = table.size();
    table.erase(std::remove_if(m_ftable, m_ltable, [&](RecPmtProp& pmt) { transformPmt(pmt); return !pmt.used; }), m_ltable);
    std::size_t fsize = table.size();

    LogDebug << isize << " -> " << fsize << " = " << isize - fsize << " PMTs are removed\n";
    return;
}

void TtCrossTalkTransformer::transformPmt(RecPmtProp& pmt) {
    if (!pmt.used || !checkPmtType(pmt)) return;
    Identifier id = TtID::id(pmt.id);
    int ch_id = TtID::getIntID(TtID::wall_id(id), TtID::lower_pmt_index(TtID::pmt(id)), TtID::strip(id));
    if (m_ch_ids.find(ch_id) == m_ch_ids.end()) return;
    pmt.used = false;
    return;
}

void TtCrossTalkTransformer::getHitsStrips() {
    m_hit_strips.clear();
    int wall_id = 0, lwr_pmt_id = 0, strip_id = 0, wall_lwr_pmt_idx = 0, ch_id_lwr_pmt = 0;
    for (RecPmtTable::const_iterator it = m_ftable; it != m_ltable; ++it) {
        Identifier id = TtID::id(it->id);
        wall_id = TtID::wall_id(id);
        lwr_pmt_id = TtID::lower_pmt_index(TtID::pmt(id));
        strip_id = TtID::strip(id);
        wall_lwr_pmt_idx = TtID::getIntID(wall_id, lwr_pmt_id, 0);
        ch_id_lwr_pmt = TtID::getIntID(wall_id, lwr_pmt_id, strip_id);

        bool found = false;
        for (std::size_t k = 0; k < m_hit_strips[wall_lwr_pmt_idx].size(); ++k) {
            if (m_hit_strips[wall_lwr_pmt_idx][k].strip_id == strip_id) {
                m_hit_strips[wall_lwr_pmt_idx][k].sum_pe += it->q;
                found = true;
                break;
            }
        }

        if (!found) {
            m_hit_strips[wall_lwr_pmt_idx].push_back({strip_id, ch_id_lwr_pmt, it->q});
        }
    }
    return;
}

void TtCrossTalkTransformer::filterCrossTalk() {
    m_ch_ids.clear();
    for (std::unordered_map<int, std::vector<HitStrip>>::iterator it = m_hit_strips.begin(); it != m_hit_strips.end(); ++it) {
        std::vector<HitStrip>& hs = it->second;
        if (hs.size() < 2) continue;
        double max_pe = 0.;
        int strip_id_max_pe = 0;
        for (const HitStrip& h : hs) {
            if (max_pe < h.sum_pe) {
                max_pe = h.sum_pe;
                strip_id_max_pe = h.strip_id;
            }
        }
        int m = - (strip_id_max_pe + 1) % 8;

        for (std::size_t k = 0; k < hs.size(); ++k) {
            if (hs[k].strip_id == strip_id_max_pe) continue;
            if (3. < hs[k].sum_pe || 0.3 < hs[k].sum_pe / max_pe) continue;

            int d_strip_id = std::abs(strip_id_max_pe - hs[k].strip_id);

            if (std::abs(d_strip_id) == 8) m_ch_ids.insert(hs[k].ic);
            else if (std::abs(d_strip_id) == 16) m_ch_ids.insert(hs[k].ic);
            else if (m != 0 && m != -1 && std::abs(d_strip_id) == 1) m_ch_ids.insert(hs[k].ic);
            else if ( (-7 < m && m < -2) && std::abs(d_strip_id) == 2) m_ch_ids.insert(hs[k].ic);
            else if ( (m == -1 || m == -2) && d_strip_id == -2) m_ch_ids.insert(hs[k].ic);
            else if ( (m == -7 || m == 0) && d_strip_id == 2) m_ch_ids.insert(hs[k].ic);
            else if( (-8 < m && m < -1) && (std::abs(d_strip_id) == 7 || std::abs(d_strip_id) == 9) ) m_ch_ids.insert(hs[k].ic);
            else if( (-8 < m && m < -1) && (std::abs(d_strip_id) == 15 || std::abs(d_strip_id) == 17) ) m_ch_ids.insert(hs[k].ic);
            else if( (-7 < m && m < -2) && (std::abs(d_strip_id) == 6 || std::abs(d_strip_id) == 10) ) m_ch_ids.insert(hs[k].ic);
            else if(m == -1 &&  (d_strip_id == 7 || d_strip_id == -9 ) ) m_ch_ids.insert(hs[k].ic);
            else if(m == -1 &&  (d_strip_id == 15 || d_strip_id == -17 ) ) m_ch_ids.insert(hs[k].ic);
            else if(m == -1 && d_strip_id == -1) m_ch_ids.insert(hs[k].ic);
            else if( (m == -1 || m == -2) && (d_strip_id == 6 || d_strip_id == -10) ) m_ch_ids.insert(hs[k].ic);
            else if(m == 0 && d_strip_id == 1) m_ch_ids.insert(hs[k].ic);
            else if(m == 0 && (d_strip_id == -7 || d_strip_id == 9) ) m_ch_ids.insert(hs[k].ic);
            else if(m == 0 && (d_strip_id == -15 || d_strip_id == 17) ) m_ch_ids.insert(hs[k].ic);
            else if( (m == 0 || m == -7) && (d_strip_id == -6 || d_strip_id == 10) ) m_ch_ids.insert(hs[k].ic);
        
        }
    }
    return;
}