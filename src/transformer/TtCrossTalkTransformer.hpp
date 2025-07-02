#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_TTCROSSTALKTRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_TTCROSSTALKTRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

#include <unordered_map>
#include <unordered_set>
#include <vector>

class TtCrossTalkTransformer : public Transformer {

    struct HitStrip {

        int strip_id;
        int ic;
        double sum_pe;

    };

public:

    TtCrossTalkTransformer(const std::string& name);

    ~TtCrossTalkTransformer() override = default;

    void operator()(RecPmtTable& table) override;

protected:

    void transform(RecPmtProp& pmt) override;
    void getHitsStrips();
    void filterCrossTalk();

    std::unordered_map<int, std::vector<HitStrip>> m_hit_strips;
    std::unordered_set<int> m_ch_ids;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_TTCROSSTALKTRANSFORMER_HPP_