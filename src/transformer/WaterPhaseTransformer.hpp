#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_WATERPHASETRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_WATERPHASETRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

#include <memory>
#include <utility>

#include <TH1D.h>

#undef __USE_WATERPHASETRANSFORMER_J24_3_0_SIMULATIONS_IMPLEMENTATION__

class WaterPhaseTransformer : public Transformer {

public:

    ~WaterPhaseTransformer() override = default;

    void transform(RecPmtTable& table) override;


#ifdef __USE_WATERPHASETRANSFORMER_J24_3_0_SIMULATIONS_IMPLEMENTATION__

    WaterPhaseTransformer(const std::string& name, double q_thold, int nb_bins, double xmin, double xmax, double lpmt_thold_min, double lpmt_thold_max, double spmt_thold_min, double spmt_thold_max, double radius, unsigned int nb_neigh_thold, double neigh_fht_diff);

protected:

    std::unique_ptr<TH1D> m_hist;

    double m_q_thold;
    double m_lpmt_thold_min, m_lpmt_thold_max;
    double m_spmt_thold_min, m_spmt_thold_max;
    double m_r2;
    double m_nb_neigh_thold;
    double m_neigh_fht_diff;
    double m_lpmt_t_i, m_lpmt_t_f;
    double m_spmt_t_i, m_spmt_t_f;
    std::size_t m_nb_spmt, m_nb_lpmt;
    std::vector<unsigned int> m_nb_neigh;
    std::vector<double> m_fht_neigh;

    void getTimes20inch(const RecPmtTable& table);
    void getTimes3inch(const RecPmtTable& table);

#else

    WaterPhaseTransformer(const std::string& name);
    WaterPhaseTransformer(const std::string& name, int nb_bins, double xmin, double xmax, double q_thold_itime, double ratio_entries_ftime, double r, unsigned int neigh_thold, double dt);

protected:

    std::unique_ptr<TH1D> m_hist;
    double m_q_thold_itime;
    double m_ratio_entries_ftime;
    double m_r2;
    unsigned int m_neigh_thold;
    double m_dt;
    std::vector<unsigned int> m_nb_neigh;

    double m_lpmt_itime, m_lpmt_ftime;

    bool getITime20inch(const RecPmtTable& table);
    bool getFTime20inch(const RecPmtTable& table);

#endif // __USE_WATERPHASETRANSFORMER_J24_3_0_SIMULATIONS_IMPLEMENTATION__

    void transformPmt(RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_WATERPHASETRANSFORMER_HPP_