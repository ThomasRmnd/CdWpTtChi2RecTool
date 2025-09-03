#include "estimator/fht/map/file/CorrectionFile.hpp"

#include "SniperKernel/SniperLog.h"

CorrectionFile::MapType CorrectionFile::s_cache;

CorrectionFile::CorrectionFile(const std::string& filename) {
    m_file = TFile::Open(filename.c_str(), "READ");
}

std::shared_ptr<CorrectionFile> CorrectionFile::open(const std::string& filename) {
    MapType::iterator it = s_cache.find(filename);
    if (it != s_cache.end()) {
        if (std::shared_ptr<CorrectionFile> existing = it->second.lock()) {
            LogDebug << "File " << filename << " already opened! Using the already existing pointer\n";
            return existing;
        }
    }
    LogDebug << "Opening file " << filename << " for the first time\n";
    std::shared_ptr<CorrectionFile> file(new CorrectionFile(filename));
    s_cache[filename] = file;
    return file;
}

bool CorrectionFile::isOpen() const {
    return (m_file && m_file->IsOpen());
}