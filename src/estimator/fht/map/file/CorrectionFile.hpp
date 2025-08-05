#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_FILE_CORRECTIONFILE_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_FILE_CORRECTIONFILE_HPP_

#include <memory>
#include <string>
#include <unordered_map>

#include <TFile.h>

/**
 * @class CorrectionFile
 * @brief Wrapper around ROOT TFile with caching and typed object access.
 *
 * This class provides a cached mechanism for opening ROOT TFiles by filename.
 * It avoids reopening the same file multiple times by maintaining a static cache.
 * Also provides typed access to objects stored in the file.
 */
class CorrectionFile {

public:
    
    /**
     * @brief Opens a ROOT file with the given filename and returns a shared pointer.
     *
     * If the file has already been opened before, returns a shared pointer to the
     * existing instance from cache. Otherwise, opens the file in READ mode.
     *
     * @param filename Path to the ROOT file.
     * 
     * @return std::shared_ptr<CorrectionFile> Shared pointer to the file handler.
     */
    static std::shared_ptr<CorrectionFile> open(const std::string& filename);

    // Checks if the file is open
    bool isOpen() const;

    /**
     * @brief Retrieves an object of type _Tp from the ROOT file by name.
     *
     * @tparam _Tp Type of the object to retrieve.
     * 
     * @param objname Name of the object in the ROOT file.
     * 
     * @return _Tp* Pointer to the requested object, or nullptr if not found.
     */
    template<typename _Tp>
    _Tp* get(const std::string& objname) const {
        return m_file->Get<_Tp>(objname.c_str());
    }

private:

    typedef std::unordered_map<std::string, std::weak_ptr<CorrectionFile>> MapType;

    TFile* m_file = nullptr; // let ROOT manage its lifetime
    static MapType s_cache;

    CorrectionFile(const std::string& filename);

};


#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_MAP_FILE_CORRECTIONFILE_HPP_