#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_CORRECTIONMAP_HPP_

#include "estimator/fht/map/file/CorrectionFile.hpp"
#include "utils/Configurable.hpp"
#include "utils/RecPmtProp.hpp"
#include "utils/TrackParams.hpp"

class CorrectionMap : public IParamsHandler, public PmtTypeChecker, public Configurable {

public:

    CorrectionMap(const std::string& name, const RecPmtType& pmt_type, const std::string& filename, const std::string& mapname);

    virtual ~CorrectionMap() = default;

    virtual bool initialize() override;

    virtual void correct(RecPmtTable::iterator first, RecPmtTable::iterator last, const double* params) = 0;

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