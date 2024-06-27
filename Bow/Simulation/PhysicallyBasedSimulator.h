#pragma once

#include <Bow/Utils/Serialization.h>
#include <Bow/Utils/FileSystem.h>
#include <Bow/Utils/Timer.h>

namespace Bow {

template <class T, int dim>
class PhysicallyBasedSimulator {
public:
    int end_frame = 120;
    T frame_dt = (T)1. / 24;
    T suggested_dt = frame_dt;
    int sub_step = 0;
    T time_elapsed = 0;
    std::string output_directory = "simulation_outputs";
    bool line_search = true;

    virtual void initialize(){};
    virtual void calculate_dt(){};
    virtual void restart_prepare(){};
    virtual void advance(T dt) = 0;
    std::function<void(void)> timestep_callback = []() {};
    virtual void dump_output(int frame_num) = 0;

    void run(int start_frame = 0)
    {
        FileSystem::create_path(output_directory);
        initialize();
        timestep_callback();
        if (!start_frame)
            dump_output(0);
        else {
            SERIALIZATION_LOAD(this->output_directory + "/restart_" + std::to_string(start_frame));
            restart_prepare();
        }
        for (int f = start_frame; f < end_frame; ++f) {
            T remaining_time = frame_dt;
            while (remaining_time) {
                calculate_dt();
                if (remaining_time <= suggested_dt) {
                    advance(remaining_time);
                    time_elapsed += remaining_time;
                    remaining_time = 0;
                }
                else if (remaining_time <= suggested_dt * 2) {
                    advance(remaining_time / 2);
                    time_elapsed += remaining_time / 2;
                    remaining_time /= 2;
                }
                else {
                    advance(suggested_dt);
                    time_elapsed += suggested_dt;
                    remaining_time -= suggested_dt;
                }
                sub_step += 1;
                timestep_callback();
                Timer::progress("Frame  " + std::to_string(f + 1), f + 1, end_frame);
                Timer::progress("Substep" + std::string(std::to_string(f + 1).length(), ' '), frame_dt - remaining_time, frame_dt);
                Timer::flush();
            }
            dump_output(f + 1);
            SERIALIZATION_SAVE(this->output_directory + "/restart_" + std::to_string(f + 1));
        }
    }
};

} // namespace Bow