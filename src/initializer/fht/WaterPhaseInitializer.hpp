#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_FHT_WATERPHASEINITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_FHT_WATERPHASEINITIALIZER_HPP_

#include "initializer/Initializer.hpp"

#include <memory>

#include <TH1D.h>

class WaterPhaseInitializer : public Initializer<FhtMethodTag> {

public:

    WaterPhaseInitializer(const std::string& name);
    WaterPhaseInitializer(const std::string& name, double dt_f2itime, double dt, double q_ratio);

    ~WaterPhaseInitializer() override = default;

    bool initiate(const vector_type& table) override;

    ParamsType getOParamsType() override;

private:

    double m_dt_f2itime;
    double m_dt;
    double m_q_ratio;
    std::unique_ptr<TH1D> m_hist;
    std::vector<RecPmtTable::const_iterator> m_it_table;

    bool getTable(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable);
    double getITime();
    vec3 getIPos();
    vec3 getFPos(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable, double itime);

};

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_FHT_WATERPHASEINITIALIZER_HPP_