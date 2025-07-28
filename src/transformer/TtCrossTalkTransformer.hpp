#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_TTCROSSTALKTRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_TTCROSSTALKTRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @class TtCrossTalkTransformer
 * 
 * @brief Derived class for transformer.
 * This transformer is used to remove cross-talk for TT PMTs
 */
class TtCrossTalkTransformer : public Transformer {

    struct HitStrip {

        int strip_id;
        int ic;
        double sum_pe;

    };

public:

    TtCrossTalkTransformer(const std::string& name);

    ~TtCrossTalkTransformer() override = default;

    /**
     * @brief Transform (remove/change) all PMTs of the target PMT type
     * 
     * @param table experiemental data vector
     */
    void transform(RecPmtTable& table) override;

protected:

    void transformPmt(RecPmtProp& pmt) override;
    void getHitsStrips();
    void filterCrossTalk();

    std::unordered_map<int, std::vector<HitStrip>> m_hit_strips;
    std::unordered_set<int> m_ch_ids;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_TTCROSSTALKTRANSFORMER_HPP_