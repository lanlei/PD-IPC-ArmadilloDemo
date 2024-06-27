#ifndef PROFILER_H
#define PROFILER_H
#include <string>
#include <chrono>
#include <Bow/Macros.h>
namespace Bow {
namespace Timer {

BOW_INLINE void flush();
BOW_INLINE void progress(std::string description, double current_t, double duration);

class ScopedTimer {
    std::chrono::time_point<std::chrono::steady_clock> m_start_time;
    std::chrono::time_point<std::chrono::steady_clock> m_end_time;
    int m_id;
    std::string m_name;
    bool m_analyze;

public:
    BOW_INLINE ScopedTimer(const std::string name, bool analyze = false);
    BOW_INLINE ~ScopedTimer();
};
}
} // namespace Bow::Timer

#ifndef BOW_STATIC_LIBRARY
#include "Timer.cpp"
#endif

#define BOW_TIMER_FLAG(name) Bow::Timer::ScopedTimer scoped_timer(name, false)
#define BOW_TIMER_ANALYZE(name) Bow::Timer::ScopedTimer scoped_timer(name, true)

#endif