#ifndef CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_

#include <string>

#include "SniperKernel/SniperJSON.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolBase.h"

namespace MyLogger {

enum Level {
    Test = 0,
    Debug = 2,
    Info = 3,
    Warn = 4,
    Error = 5,
    Fatal = 6,
};

#define MYSNIPERLOG(Flag, ObjName) SniperLog::Logger(Flag, scope(), ObjName, __func__)
#define MYLOG(Flag, ObjName) (logLevel() > Flag ? SniperLog::Logger::Silencer : MYSNIPERLOG(Flag, ObjName))

SniperLog::Logger Log(Level level, const char* objname) {
    return MYLOG(level, objname);
}

} // namespace MyLogger

class Configurable : public ToolBase {

public:

    Configurable() : ToolBase{""} {}
    Configurable(const std::string& name) : ToolBase{name} {}

    virtual ~Configurable() = default;

    virtual void configure(const SniperJSON&) {}

protected:

    template<typename _Tp>
    bool setConfigValue(_Tp& member, const std::string& membername, const SniperJSON& config) {
        std::string fullkey = m_name + "__" + membername;
        SniperJSON::map_iterator it = config.find(fullkey);
        if (it == config.map_end()) return false;
        member = it->second.get<_Tp>();
        LogDebug << membername << " = " << member << '\n';
        MyLogger::Log(MyLogger::Debug, m_name) << membername << " = " << member << '\n';
        return true;
    }

};

#endif // CDWPTTCHI2RECTOOL_UTILS_CONFIGURABLE_HPP_