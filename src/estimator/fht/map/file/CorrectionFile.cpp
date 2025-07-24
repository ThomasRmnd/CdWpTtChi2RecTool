#include "estimator/fht/map/file/CorrectionFile.hpp"

#include <TROOT.h>

#include "SniperKernel/SniperLog.h"

CorrectionFile::MapType CorrectionFile::s_cache;

CorrectionFile::CorrectionFile(const std::string& filename) {
    /* TFile* f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        LogError << "Failed to open file: " << filename;
        m_file.reset();
        return;
    }
    gROOT->GetListOfFiles()->Remove(f);
    m_file.reset(f); */
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
