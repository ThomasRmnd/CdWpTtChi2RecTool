#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_WATERPHASETRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_WATERPHASETRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

#include <memory>
#include <utility>

#include <TH1D.h>

class WaterPhaseTransformer : public Transformer {

public:

    WaterPhaseTransformer(const std::string& name, int nb_bins, double xmin, double xmax, double q_thold_itime, double ratio_entries_ftime, double r, unsigned int neigh_thold, double dt);

    ~WaterPhaseTransformer() override = default;

    void transform(RecPmtTable& table) override;

protected:

    std::unique_ptr<TH1D> m_hist;
    double m_q_thold_itime;
    double m_ratio_entries_ftime;
    double m_r2;
    unsigned int m_neigh_thold;
    double m_dt;
    std::vector<unsigned int> m_nb_neigh;

    double m_lpmt_itime, m_lpmt_ftime;

    void transformPmt(RecPmtProp& pmt) override;

    bool getITime20inch(const RecPmtTable& table);
    bool getFTime20inch(const RecPmtTable& table);

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_WATERPHASETRANSFORMER_HPP_