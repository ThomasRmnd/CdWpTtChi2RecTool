#include "estimator/fht/map/file/CorrectionFile.hpp"

CorrectionFile::MapType CorrectionFile::s_cache;

CorrectionFile::CorrectionFile(const std::string& filename) {
    m_file.reset(TFile::Open(filename.c_str(), "READ"));
}

std::shared_ptr<CorrectionFile> CorrectionFile::open(const std::string& filename) {
    MapType::iterator it = s_cache.find(filename);
    if (it != s_cache.end()) {
        if (std::shared_ptr<CorrectionFile> existing = it->second.lock()) {
            return existing;
        }
    }
    std::shared_ptr<CorrectionFile> file = std::make_shared<CorrectionFile>(filename);
    s_cache[filename] = file;
    return file;
}

bool CorrectionFile::isOpen() const {
    return (m_file && m_file->IsOpen());
}

template<typename _Tp>
_Tp* CorrectionFile::get(const std::string& objname) const {
    return m_file->Get<_Tp>(objname.c_str());
}
