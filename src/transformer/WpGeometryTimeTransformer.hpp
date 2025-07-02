#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_WPGEOMETRYTIMETRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_WPGEOMETRYTIMETRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

#include <memory>

#include <TH1D.h>

struct ClusteredPmt {

    RecPmtTable::const_iterator it_pmt;
    int cluster_id;

};

class WpGeometryTimeTransformer : public Transformer {

public:

    WpGeometryTimeTransformer(const std::string& name);
    WpGeometryTimeTransformer(const std::string& name, double time_window_early, double time_window_late, double thold_sep, double thold_q_ratio, double radius_time, double radius_arclength);

    ~WpGeometryTimeTransformer() override = default;

    void operator()(RecPmtTable& table) override;

protected:

    double m_time_window_early;
    double m_time_window_late;
    
    double m_thold_sep;
    double m_thold_q_ratio;

    unsigned int m_thold_nb_neighbors;
    double m_radius_time;
    double m_radius_arclength;

    double m_time_early, m_time_late;
    vec3 m_pos_early, m_pos_late;
    vec3 m_plan_z;
    double m_max_range;
    double m_max_q;
    std::vector<ClusteredPmt> m_clustered_table;

    //
    std::vector<unsigned int> m_nb_neigh;
    double m_dist2_neigh;
    unsigned int m_neigh_thold;
    double m_pmt_thold;
    std::unique_ptr<TH1D> m_h_fht;
    double m_first_time_late;
    double m_width_time_window;
    //

    // //
    double m_qthold_iftime;
    // //

    double m_dist_range, m_time_range;

    bool getEarlyLateTime();
    bool getEarlyLatePosition();
    bool getMaxRange();
    bool fillClusteredTable();
    void dbscan(unsigned int thold_nb_neighbors);
    void expandCluster(std::vector<std::size_t>& neighbors, int cluster_id, unsigned int thold_nb_neighbors);
    std::vector<std::size_t> regionQuery(const std::size_t idx);
    void checkClusters();
    int getClusterWithMaxPmts();
    bool getMaxChargeInCluster(int cluster_id);
    void getMeanPosAndTimeInCluster(int cluster_id);
    double arcLengthRatio(const RecPmtProp& pmt, double range);

    void transform(RecPmtProp& pmt) override;

};


#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_WPGEOMETRYTIMETRANSFORMER_HPP_