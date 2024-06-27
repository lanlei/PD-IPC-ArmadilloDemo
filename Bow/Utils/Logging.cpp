#include "Logging.h"
#include <iostream>
#include <string>
#include <fstream>
#include <mutex>
#include <unordered_map>
#include <chrono>
#include <vector>
#include <memory>

/*
 * https://gist.github.com/kevinkreiser/39f2e39273c625d96790
 */

/*
ATTRIBUTE_RESET = "\033[0m";
BLACK = "\033[22;30m";
RED = "\033[22;31m";
GREEN = "\033[22;32m";
BROWN = "\033[22;33m";
BLUE = "\033[22;34m";
MAGENTA = "\033[22;35m";
CYAN = "\033[22;36m";
GREY = "\033[22;37m";
DARKGREY = "\033[01;30m";
LIGHTRED = "\033[01;31m";
LIGHTGREEN = "\033[01;32m";
YELLOW = "\033[01;33m";
LIGHTBLUE = "\033[01;34m";
LIGHTMAGENTA = "\033[01;35m";
LIGHTCYAN = "\033[01;36m";
WHITE = "\033[01;37m";
*/

namespace Bow {
namespace Logging {
namespace internal {

const std::unordered_map<LogLevel, std::string> log_flag{
    { Timing, "\x1b[36;1m[Timing]\x1b[0m " },
    { Fatal, "\x1b[31;1m[FATAL]\x1b[0m " },
    { Error, "\x1b[31;1m[ERROR]\x1b[0m " },
    { Warning, "\x1b[33;1m[WARN]\x1b[0m " },
    { Info, "\x1b[32;1m[INFO]\x1b[0m " },
    { Debug, "\x1b[34;1m[DEBUG]\x1b[0m " }
};

const std::unordered_map<LogLevel, std::string> log_flag_uncolored{
    { Timing, "[Timing] " },
    { Fatal, "[FATAL] " },
    { Error, "[ERROR] " },
    { Warning, "[WARN] " },
    { Info, "[INFO] " },
    { Debug, "[DEBUG] " }
};

using LogConfig = std::unordered_map<std::string, std::string>;
class Logger {
public:
    Logger() = delete;
    Logger(const LogLevel level)
        : m_level(level) {}
    void set_level(const LogLevel level) { m_level = level; }
    virtual ~Logger(){};
    virtual void log(const std::string&, const LogLevel level) = 0;

protected:
    virtual void log(const std::string&) = 0;
    std::mutex m_lock;
    LogLevel m_level;
};

class StdOutLogger : public Logger {
public:
    StdOutLogger() = delete;
    StdOutLogger(const LogLevel level)
        : Logger(level) {}
    inline virtual void log(const std::string& message, const LogLevel level)
    {
        if (level < Logger::m_level)
            return;
        std::string output;
        output.reserve(message.length() + 64);
        output.append(log_flag.at(level));
        output.append(message);
        output.push_back('\n');
        log(output);
    }

protected:
    inline virtual void log(const std::string& message)
    {
        std::cout << message;
        std::cout.flush();
    }
};

class FileLogger : public Logger {
public:
    FileLogger() = delete;
    FileLogger(const LogLevel level, const std::string& filename)
        : Logger(level)
    {
        m_file_name = filename;
        reopen();
    }
    ~FileLogger()
    {
        try {
            m_file.close();
        }
        catch (...) {
        }
    }
    inline virtual void log(const std::string& message, const LogLevel level)
    {
        if (level < m_level)
            return;
        std::string output;
        output.reserve(message.length() + 64);
        //   output.append(timestamp());
        output.append(log_flag_uncolored.at(level));
        output.append(message);
        output.push_back('\n');
        log(output);
    }

protected:
    inline virtual void log(const std::string& message)
    {
        Logger::m_lock.lock();
        m_file << message;
        m_file.flush();
        Logger::m_lock.unlock();
        reopen();
    }

    inline void reopen()
    {
        //check if it should be closed and reopened
        auto now = std::chrono::system_clock::now();
        Logger::m_lock.lock();
        if (now - m_last_reopen > m_reopen_interval) {
            m_last_reopen = now;
            try {
                m_file.close();
            }
            catch (...) {
            }
            try {
                m_file.open(m_file_name, std::ofstream::out | std::ofstream::app);
                m_last_reopen = std::chrono::system_clock::now();
            }
            catch (std::exception& e) {
                try {
                    m_file.close();
                }
                catch (...) {
                }
                throw e;
            }
        }
        Logger::m_lock.unlock();
    }
    std::string m_file_name;
    std::ofstream m_file;
    std::chrono::seconds m_reopen_interval = std::chrono::seconds(300);
    std::chrono::system_clock::time_point m_last_reopen;
};

static std::vector<std::unique_ptr<Logger>> logger_pool;

} // namespace internal

BOW_INLINE void set_level(const LogLevel level)
{
    if (internal::logger_pool.size() == size_t(0)) {
        internal::logger_pool.emplace_back(std::make_unique<internal::StdOutLogger>(Info));
    }
    internal::logger_pool[0]->set_level(level);
}
BOW_INLINE void log_message(const std::string message, const LogLevel level)
{
    if (internal::logger_pool.size() == size_t(0)) {
        internal::logger_pool.emplace_back(std::make_unique<internal::StdOutLogger>(Info));
    }
    for (auto& logger_ptr : internal::logger_pool) {
        logger_ptr->log(message, level);
    }
    if (level == Fatal)
        throw std::runtime_error(message.c_str());
}
BOW_INLINE void new_logger(const std::string filename, const LogLevel level, bool reset_file)
{
    if (internal::logger_pool.size() == size_t(0)) {
        internal::logger_pool.emplace_back(std::make_unique<internal::StdOutLogger>(Info));
    }
    if (reset_file) {
        std::ofstream file(filename);
        file.close();
    }

    internal::logger_pool.emplace_back(std::make_unique<internal::FileLogger>(level, filename));
}
}
} // namespace Bow::Logging