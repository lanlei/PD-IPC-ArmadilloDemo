#ifndef LOGGING_H
#define LOGGING_H
#include <string>
#include <Bow/Macros.h>
#include <vector>
#include <sstream>

namespace Bow {
namespace Logging {
enum LogLevel {
    Debug = 0,
    Info = 100,
    Timing = 200,
    Warning = 500,
    Error = 1000,
    Fatal = 9000,
    None = 9001
};
/*
 * Only change the level of the StdOutLogger.
 */
BOW_INLINE void set_level(const LogLevel level);
/*
 * Log a message. The message will be logged into all loggers in the logger pool.
 * By default, a StdOutLogger will be automatically created. The message will be always shown in the terminal.
 */
BOW_INLINE void log_message(const std::string message, const LogLevel level);
/*
 * Push a new file logger with a given level into the logger pool.
 */
BOW_INLINE void new_logger(const std::string filename, const LogLevel level = Info, bool reset_file = false);

inline void collect_message(std::stringstream& ss)
{
    return;
}

template <typename U, typename... Ts>
inline void collect_message(std::stringstream& ss, const U& first_obj, const Ts&... objs)
{
    ss << first_obj;
    collect_message(ss, objs...);
}

template <typename... Ts>
inline void log(const LogLevel level, const Ts&... objs)
{
    std::stringstream ss;
    collect_message(ss, objs...);
    log_message(ss.str(), level);
}

template <typename... Ts>
inline void debug(const Ts&... objs)
{
    log(Debug, objs...);
}
template <typename... Ts>
inline void info(const Ts&... objs)
{
    log(Info, objs...);
}
template <typename... Ts>
inline void timing(const Ts&... objs)
{
    log(Timing, objs...);
}
template <typename... Ts>
inline void warn(const Ts&... objs)
{
    log(Warning, objs...);
}
template <typename... Ts>
inline void error(const Ts&... objs)
{
    log(Error, objs...);
}
template <typename... Ts>
inline void fatal(const Ts&... objs)
{
    log(Fatal, objs...);
}
} // namespace Logging

#define assert_info(x, info)                   \
    {                                          \
        bool ___ret___ = static_cast<bool>(x); \
        if (!___ret___) {                      \
            Logging::error(info);              \
            exit(0);                           \
        }                                      \
    }

#define BOW_ASSERT(x) BOW_ASSERT_INFO((x), #x)
#define BOW_ASSERT_INFO assert_info
#define BOW_NOT_IMPLEMENTED BOW_ASSERT_INFO(false, std::string(__func__) + " not implemented.")

} // namespace Bow

#ifndef BOW_STATIC_LIBRARY
#include "Logging.cpp"
#endif

#endif