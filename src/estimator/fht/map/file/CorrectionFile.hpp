#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_FILE_CORRECTIONFILE_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_FILE_CORRECTIONFILE_HPP_

#include <memory>
#include <string>
#include <unordered_map>

#include <TFile.h>

class CorrectionFile {

public:
    
    static std::shared_ptr<CorrectionFile> open(const std::string& filename);

    bool isOpen() const;

    template<typename _Tp>
    _Tp* get(const std::string& objname) const {
        return m_file->Get<_Tp>(objname.c_str());
    }

private:

    typedef std::unordered_map<std::string, std::weak_ptr<CorrectionFile>> MapType;

    TFile* m_file = nullptr;
    static MapType s_cache;

    CorrectionFile(const std::string& filename);

};


#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_FILE_CORRECTIONFILE_HPP_