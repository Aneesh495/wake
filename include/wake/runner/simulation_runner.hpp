#pragma once

#include "wake/scenarios/scenario.hpp"
#include "wake/runner/thread_pool.hpp"
#include <memory>
#include <atomic>
#include <functional>
#include <string>

namespace wake {

struct RunnerConfig {
    uint64_t max_steps{10000};
    uint64_t export_interval{200};
    uint64_t render_interval{50};
    bool live_terminal_render{true};
    std::string export_dir{"./output"};
    double target_fps{60.0};
};

class SimulationRunner {
public:
    explicit SimulationRunner(
        std::shared_ptr<Scenario> scenario,
        const RunnerConfig& config = RunnerConfig()
    );

    ~SimulationRunner();

    // Blocking synchronous execution loop until scenario completion or signal
    void run();

    // Asynchronous non-blocking background execution
    void start_async();
    void pause();
    void resume();
    void stop();
    void wait();

    bool is_running() const noexcept { return is_running_; }
    bool is_paused() const noexcept { return is_paused_; }
    uint64_t current_step() const;
    double average_mlups() const noexcept { return avg_mlups_; }

    // User callback hooks
    void set_step_callback(std::function<void(uint64_t step, double mlups)> callback) {
        step_callback_ = std::move(callback);
    }

private:
    void execution_loop();

    std::shared_ptr<Scenario> scenario_;
    RunnerConfig config_;

    std::atomic<bool> is_running_{false};
    std::atomic<bool> is_paused_{false};
    std::atomic<bool> should_stop_{false};

    std::thread worker_thread_;
    double avg_mlups_{0.0};
    uint64_t step_count_{0};

    std::function<void(uint64_t step, double mlups)> step_callback_;
};

} // namespace wake
