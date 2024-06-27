#include "Timer.h"
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "Logging.h"
#include <algorithm>

namespace Bow {
	namespace Timer {
		namespace internal {
			static std::map<std::pair<std::string, int>, int> name_scope;
			static std::vector<std::pair<std::string, int>> scope_name(1, std::make_pair("Global", -1));
			static std::vector<std::chrono::duration<double>> scope_duration(1, std::chrono::duration<double>(0));
			static std::vector<std::chrono::duration<double>> global_duration(1, std::chrono::duration<double>(0));
			static std::vector<std::vector<int>> scope_edges(1, std::vector<int>());
			static std::vector<int> scope_stack(1, 0);

			static std::map<std::string, double> analyze_min;
			static std::map<std::string, double> analyze_max;
			static std::map<std::string, double> analyze_total;
			static std::map<std::string, int> analyze_number;
			inline void analyze_timing(std::string name, double value)
			{
				++analyze_number[name];
				if (analyze_number[name] == 1) {
					analyze_min[name] = value;
					analyze_max[name] = value;
				}
				analyze_min[name] = std::min(analyze_min[name], value);
				analyze_max[name] = std::max(analyze_max[name], value);
				analyze_total[name] += value;
				printf("%s :\n", name.c_str());
				printf("Min %.20f\n", analyze_min[name]);
				printf("Max %.20f\n", analyze_max[name]);
				printf("Average %.20f\n", analyze_total[name] / analyze_number[name]);
				printf("Total %.20f\n\n", analyze_total[name]);
			}

			class GlobalTimer {
			public:
				std::chrono::time_point<std::chrono::steady_clock> current_time;

				GlobalTimer()
				{
					current_time = std::chrono::steady_clock::now();
				}
			};
			static internal::GlobalTimer global_timer;

			inline std::string duration2string(const std::chrono::duration<double>& elapsed_seconds)
			{
				using namespace std::chrono;
				std::string str;
				auto dur = elapsed_seconds;

				int num_days = int(dur / hours(24));
				if (num_days) str += std::to_string(num_days) + "d ";
				dur -= num_days * hours(24);

				int num_hours = int(dur / hours(1));
				if (num_hours) str += std::to_string(num_hours) + "h ";
				dur -= num_hours * hours(1);

				int num_minutes = int(dur / minutes(1));
				if (num_minutes) str += std::to_string(num_minutes) + "m ";
				dur -= num_minutes * minutes(1);

				str += std::to_string(dur.count()) + "s";
				return str;
			}

			inline void traverseScopes(int id, int depth, std::stringstream& ss)
			{
				auto printScope = [&]() {
					std::string scope_str = duration2string(scope_duration[id]);
					std::string global_str = duration2string(global_duration[id]);
					double scope_percent = scope_duration[id].count() * 100 / scope_duration[0].count();
					double global_percent = global_duration[id].count() * 100 / global_duration[0].count();
					ss << std::fixed << std::setprecision(1);
					ss << " " << scope_name[id].first << " : " << scope_str << " (" << scope_percent << "%)   " << global_str << " (" << global_percent << "%)";
					Logging::timing(ss.str());
					ss.str("");
				};
				if (scope_edges[id].empty()) {
					for (int i = 0; i < depth; ++i) ss << "|||";
					printScope();
				}
				else {
					for (int i = 0; i < depth; ++i) ss << "|||";
					ss << "<<";
					printScope();
					for (auto s : scope_edges[id])
						traverseScopes(s, depth + 1, ss);

					// other scope
					std::chrono::duration<double> other_scope = scope_duration[id];
					std::chrono::duration<double> other_global = global_duration[id];
					for (auto s : scope_edges[id]) {
						other_scope -= scope_duration[s];
						other_global -= global_duration[s];
					}
					double other_scope_percent = other_scope.count() * 100 / scope_duration[0].count();
					double other_global_percent = other_global.count() * 100 / global_duration[0].count();
					if (other_global_percent >= 0.05) {
						for (int i = 0; i <= depth; ++i) ss << "|||";
						std::string scope_str = duration2string(other_scope);
						std::string global_str = duration2string(other_global);
						ss << " Uncounted : " << scope_str << " (" << other_scope_percent << "%)   " << global_str << " (" << other_global_percent << "%)";
						Logging::timing(ss.str());
						ss.str("");
					}

					for (int i = 0; i < depth; ++i) ss << "|||";
					ss << ">>";
					Logging::timing(ss.str());
					ss.str("");
				}
			}
		} // namespace internal

		BOW_INLINE ScopedTimer::ScopedTimer(const std::string name, bool analyze)
			: m_name(name)
			, m_analyze(analyze)
		{
			m_start_time = std::chrono::steady_clock::now();
			auto name_parent = std::make_pair(name, internal::scope_stack.back());
			if (internal::name_scope.find(name_parent) != internal::name_scope.end()) {
				m_id = internal::name_scope[name_parent];
			}
			else {
				m_id = internal::scope_name.size();
				internal::name_scope[name_parent] = m_id;
				internal::scope_name.push_back(name_parent);
				internal::scope_duration.emplace_back(0);
				internal::global_duration.emplace_back(0);
				internal::scope_edges.emplace_back();
				internal::scope_edges[internal::scope_stack.back()].push_back(m_id);
			}
			internal::scope_stack.push_back(m_id);
		}

		BOW_INLINE ScopedTimer::~ScopedTimer()
		{
			m_end_time = std::chrono::steady_clock::now();
			internal::scope_duration[m_id] += m_end_time - m_start_time;
			internal::scope_stack.pop_back();
			if (m_analyze) {
				internal::analyze_timing(m_name + std::string(" Timing"), (m_end_time - m_start_time).count() / 1e9);
			}
		}

		BOW_INLINE void progress(std::string description, double current_t, double duration)
		{
			int length = 50;
			std::stringstream ss;
			ss << description << " ";
			ss << "[";
			double simu_percent = (double)current_t / duration;
			int simu_int = (int)(simu_percent * length);
			for (int i = 0; i < simu_int; ++i) ss << "#";
			for (int i = simu_int; i < length; ++i) ss << " ";
			ss << "]";
			ss << std::fixed << std::setprecision(2) << std::setw(8);
			ss << current_t / duration * 100 << "%";
			Logging::info(ss.str());
		}

		BOW_INLINE void flush()
		{
			std::chrono::time_point<std::chrono::steady_clock> last_time = internal::global_timer.current_time;
			std::chrono::time_point<std::chrono::steady_clock> current_time = std::chrono::steady_clock::now();
			internal::global_timer.current_time = current_time;
			internal::scope_duration[0] = current_time - last_time;
			for (size_t i = 0; i < internal::scope_duration.size(); ++i)
				internal::global_duration[i] += internal::scope_duration[i];
			std::stringstream ss;
			internal::traverseScopes(0, 0, ss);
			// std::cout << ss.str();
			for (size_t i = 0; i < internal::scope_duration.size(); ++i)
				internal::scope_duration[i] = std::chrono::duration<double>(0);
			puts("");
		}
	}
} // namespace Bow::Timer