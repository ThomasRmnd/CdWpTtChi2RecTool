#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERMAXCHARGEINITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERMAXCHARGEINITIALIZER_HPP_

#include "initializer/Initializer.hpp"

#include <memory>

#include <TH1D.h>

class ClusterMaxChargeInitializer : public Initializer<FhtMethodTag> {

public:

    ClusterMaxChargeInitializer(const std::string& name);
    ClusterMaxChargeInitializer(const std::string& name, double pmt_cnt_thold, double shift, double range, double pe_thold, double radius, double lwr_q_thold, double dir_corr_factor, double radius_end);

    ~ClusterMaxChargeInitializer() override = default;

    bool initiate(const RecPmtTable& table) override;

    ParamsType getOParamsType() override;

private:

    double m_pmt_cnt_thold;
    double m_shift;
    double m_range;
    double m_pe_thold;
    double m_radius;
    double m_lwr_q_thold;
    double m_dir_corr_factor;
    double m_radius_end;
    std::unique_ptr<TH1D> m_hist;

    std::vector<RecPmtTable::const_iterator> m_it_table;

    bool getFhtCharge(RecPmtTable::const_iterator it, float& q);
    bool getTable(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable);
    double getITime();
    vec3 getIPos();
    vec3 getChargeCentroid(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable);
    vec3 getFPos(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable, const vec3& pos_i);

};

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERMAXCHARGEINITIALIZER_HPP_