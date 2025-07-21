#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP_HPP_

#include <TFile.h>

#include "SniperKernel/ToolBase.h"

#include "estimator/fht/map/file/CorrectionFile.hpp"
#include "utils/RecPmtProp.hpp"
#include "utils/TrackParams.hpp"
#include "utils/vec3.hpp"

template<ParamsType _Pt>
class CorrectionMap : public ToolBase, public PmtTypeChecker, public TrackSetter<_Pt> {

public:

    CorrectionMap(const std::string& name) :
        ToolBase(name),
        PmtTypeChecker(RecPmtType::PMT_UNKNOWN)
    {}

    CorrectionMap(const std::string& name, const RecPmtType& pmt_type, const std::string& filename, const std::string& mapname) :
        ToolBase(name),
        PmtTypeChecker(pmt_type),
        m_filename(filename),
        m_mapname(mapname)
    {}

    virtual ~CorrectionMap() = default;

    virtual bool initialize() override {
        if (!openCorrFile()) return false;
        if (!openCorrProfile()) return false;
        return true;
    }

    virtual void operator()(RecPmtTable& table) = 0;

protected:

    std::shared_ptr<CorrectionFile> m_file = nullptr;
    std::string m_filename;
    std::string m_mapname;

    bool openCorrFile() {
        m_file = CorrectionFile::open(m_filename);
        if (!m_file->isOpen()) {
            LogError << "Cannot open correction file " << m_filename << '\n';
            return false;
        }
        return true;
    }

    virtual bool openCorrProfile() = 0;

    virtual double correction(const RecPmtProp& pmt) = 0;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP_HPP_