#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_HYBRIDLSWATERCORRECTIONMAP3D_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_HYBRIDLSWATERCORRECTIONMAP3D_HPP_

#include <TProfile3D.h>

#include "estimator/fht/map/CorrectionMap.hpp"
#include "estimator/fht/map/CorrParam.hpp"

class HybridLSWaterCorrectionMap3d : public CorrectionMap, public TrackSetter<SingleCdParamsTag> {

public:

    HybridLSWaterCorrectionMap3d(const std::string& name, const RecPmtType& pmt_type, const std::string& filename, const std::string& mapname, const std::shared_ptr<CorrParam>& param_x, const std::shared_ptr<CorrParam>& param_y, const std::shared_ptr<CorrParam>& param_z);

    ~HybridLSWaterCorrectionMap3d() override = default;

    ParamsType getIParamsType() override;

    bool initialize() override;

    void correct(RecPmtTable::iterator first, RecPmtTable::iterator last, const double* params) override;

private:

    std::shared_ptr<CorrParam> m_param_x;
    std::shared_ptr<CorrParam> m_param_y;
    std::shared_ptr<CorrParam> m_param_z;
    TProfile3D* m_prof3d;

    int m_nbins_x, m_nbins_y, m_nbins_z;
    double m_fbin_ctr_x, m_fbin_ctr_y, m_fbin_ctr_z;
    double m_lbin_ctr_x, m_lbin_ctr_y, m_lbin_ctr_z;

    bool openCorrProfile() override;

    int getBinX(double value);
    int getBinY(double value);
    int getBinZ(double value);

    double correction(const RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_HYBRIDLSWATERCORRECTIONMAP3D_HPP_