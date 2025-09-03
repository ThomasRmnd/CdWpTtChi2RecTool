#include "estimator/fht/map/CorrectionMap.hpp"

CorrectionMap::CorrectionMap(const std::string& name, const RecPmtType& pmt_type, const std::string& filename, const std::string& mapname) :
    Configurable{name},
    PmtTypeChecker{pmt_type},
    m_filename{filename},
    m_mapname{mapname}
{}

bool CorrectionMap::initialize() {
    if (!openCorrFile()) return false;
    if (!openCorrProfile()) return false;
    return true;
}

bool CorrectionMap::openCorrFile() {
    m_file = CorrectionFile::open(m_filename);
    if (!m_file->isOpen()) {
        LogError << "Cannot open correction file " << m_filename << '\n';
        return false;
    }
    return true;
}