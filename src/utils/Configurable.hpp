#ifndef CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_

#include "SniperKernel/ToolBase.h"

#include "SniperKernel/SniperJSON.h"

class Configurable : public ToolBase {

public:

    using ToolBase::ToolBase;

    virtual ~Configurable() = default;

    virtual void configure(const SniperJSON& config) {}

protected:

    template<typename _Tp>
    _Tp getConfigValue(const std::string& membername, const SniperJSON& config) const {
        std::string fullkey = m_name + "__" + membername;
        SniperJSON::map_iterator it = config.find(fullkey);
        if (it == config.map_end()) {
            LogError << "Cannot retrieve the member value at " << fullkey << '\n';
            return _Tp{};
        }
        return it->second.get<_Tp>();
    }

};

#endif // CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_