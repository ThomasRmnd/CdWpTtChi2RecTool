#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP_HPP_

#include <TFile.h>

#include "SniperKernel/ToolBase.h"

#include "utils/RecPmtProp.hpp"
#include "utils/TrackParams.hpp"
#include "utils/vec3.hpp"

template<typename _Pt>
class CorrectionMap : public ToolBase, public PmtTypeChecker, public TrackSetter<_Pt> {

public:

    using ToolBase::ToolBase;

    CorrectionMap(const std::string& name, const PmtType& pmt_type, const std::string& corrfilename, const std::string& corrmap_name) :
        ToolBase(name),
        PmtTypeChecker(pmt_type),
        m_corrfilename(corrfilename),
        m_corrmap_name(corrmap_name),
        m_corrfile(nullptr),
        m_is_initialized(false)
    {};

    virtual ~CorrectionMap() = default;

    virtual bool initialize() override {
        if (m_is_initialized) return true;
        if (!openCorrFile()) return false;
        if (!openCorrProfile()) return false;
        return (m_is_initialized = true);
    };

    virtual bool finalize() override {
        if (m_corrfile->IsOpen()) m_corrfile->Close();
        return true;
    };

    virtual void operator()(RecPmtTable& table) = 0;

protected:

    std::string m_corrfilename;
    std::string m_corrmap_name;
    TFile* m_corrfile;
    bool m_is_initialized;

    bool openCorrFile() {
        m_corrfile = TFile::Open(m_corrfilename.c_str());
        if (!m_corrfile || !m_corrfile->IsOpen()) {
            LogError << "Cannot open correction file " << m_corrfilename << '\n';
            return false;
        }
        LogDebug << "Correction file " << m_corrfilename << " is opened\n";
        return true;
    };

    virtual bool openCorrProfile() = 0;

    virtual double correction(const RecPmtProp& pmt) = 0;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP_HPP_