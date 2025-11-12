#include "estimator/fht/map/CorrectionMap.hpp"

#include "DataPathHelper/Path.hh"

CorrectionMap::CorrectionMap(const std::string& name, const RecPmtType& pmt_type, const std::string& filename, const std::string& mapname) :
    PmtTypeChecker{pmt_type},
    Configurable{name},
    m_filename{filename},
    m_mapname{mapname}
{}

void CorrectionMap::configure(const SniperJSON& config) {
    setConfigValue(m_filename, "Filename", config);
    setConfigValue(m_mapname, "Mapname", config);
}

bool CorrectionMap::initialize() {
    if (!openCorrFile()) return false;
    if (!openCorrProfile()) return false;
    return true;
}

bool CorrectionMap::openCorrFile() {
    std::string filename = JUNO::Path::resolve(m_filename);
    m_file = CorrectionFile::open(filename);
    if (!m_file->isOpen()) {
        LogError << "Failed to open correction file: " << m_filename << "\n";
        return false;
    }
    return true;
}