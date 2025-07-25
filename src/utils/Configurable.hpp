#ifndef CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_

#include <string>

#include "SniperKernel/SniperJSON.h"
#include "SniperKernel/SniperLog.h"

class Configurable {

public:

    Configurable() : c_name{} {}
    Configurable(const std::string& name) : c_name{name} {}

    virtual ~Configurable() = default;

    virtual void configure(const SniperJSON&) {}

protected:

    const std::string c_name;

    template<typename _Tp>
    _Tp getConfigValue(const std::string& membername, const SniperJSON& config) {
        std::string fullkey = c_name + "__" + membername;
        SniperJSON::map_iterator it = config.find(fullkey);
        if (it == config.map_end()) {
            LogError << "Cannot retrieve the member value " << membername << " in " << c_name << '\n';
            return _Tp{};
        }
        _Tp value = it->second.get<_Tp>();
        LogDebug << membername << " = " << value << '\n';
        return value;
    }

};

#endif // CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_