#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERMAXCHARGEINITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_FHT_CLUSTERMAXCHARGEINITIALIZER_HPP_

#include <memory>

#include <TH1D.h>

#include "initializer/Initializer.hpp"

/**
 * @class ClusterMaxChargeInitializer
 * 
 * @brief Derived class to calculate track parameters initial guess for FHT method
 * This initial guess is designed for LS phase single track 
 */
class ClusterMaxChargeInitializer : public Initializer<FhtMethodTag> {

public:

    ClusterMaxChargeInitializer(const std::string& name, double pmt_cnt_thold, double shift, double range, double pe_thold, double ipos_radius, double q_ratio, double dir_corr_factor, double fpos_radius);

    ~ClusterMaxChargeInitializer() override = default;

    void configure(const SniperJSON& config);

    ParamsType getOParamsType() override;

    /**
     * @brief Calculate the track parameters initial guess
     * 
     * @param data experimental data vector
     * 
     * @return Boolean whether the calculation is successfull or not
     */
    bool initiate(const RecPmtTable& table) override;

private:

    double m_pmt_cnt_thold;
    double m_shift;
    double m_range;
    double m_pe_thold;
    double m_ipos_radius;
    double m_q_ratio;
    double m_dir_corr_factor;
    double m_fpos_radius;
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