#include "wake/runner/thread_pool.hpp"
#include "wake/runner/simulation_runner.hpp"
#include <iostream>
#include <chrono>
#include <filesystem>
#include <cmath>

namespace wake {

// ThreadPool Implementation
ThreadPool::ThreadPool(size_t threads) {
    if (threads == 0) threads = 1;

    for (size_t i = 0; i < threads; ++i) {
        workers_.emplace_back([this]() {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(queue_mutex_);
                    cv_task_.wait(lock, [this]() {
                        return stop_ || !tasks_.empty();
                    });

                    if (stop_ && tasks_.empty()) {
                        return;
                    }

                    task = std::move(tasks_.front());
                    tasks_.pop();
                    active_tasks_++;
                }

                task();

                {
                    std::unique_lock<std::mutex> lock(queue_mutex_);
                    active_tasks_--;
                    if (tasks_.empty() && active_tasks_ == 0) {
                        cv_finished_.notify_all();
                    }
                }
            }
        });
    }
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        stop_ = true;
    }
    cv_task_.notify_all();

    for (std::thread& worker : workers_) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

void ThreadPool::wait_idle() {
    std::unique_lock<std::mutex> lock(queue_mutex_);
    cv_finished_.wait(lock, [this]() {
        return tasks_.empty() && active_tasks_ == 0;
    });
}

void ThreadPool::parallel_for(
    int start,
    int end,
    const std::function<void(int start_chunk, int end_chunk)>& chunk_fn
) {
    int total_work = end - start;
    if (total_work <= 0) return;

    size_t num_threads = workers_.size();
    if (num_threads == 0) num_threads = 1;

    int chunk_size = (total_work + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);
    std::vector<std::future<void>> futures;

    for (size_t t = 0; t < num_threads; ++t) {
        int c_start = start + static_cast<int>(t) * chunk_size;
        int c_end = std::min(end, c_start + chunk_size);
        if (c_start < c_end) {
            futures.push_back(enqueue([c_start, c_end, &chunk_fn]() {
                chunk_fn(c_start, c_end);
            }));
        }
    }

    for (auto& f : futures) {
        f.get();
    }
}

// SimulationRunner Implementation
SimulationRunner::SimulationRunner(
    std::shared_ptr<Scenario> scenario,
    const RunnerConfig& config
) : scenario_(std::move(scenario)),
    config_(config) {}

SimulationRunner::~SimulationRunner() {
    stop();
    if (worker_thread_.joinable()) {
        worker_thread_.join();
    }
}

void SimulationRunner::run() {
    if (!scenario_) {
        throw std::runtime_error("Cannot run null scenario in SimulationRunner");
    }

    is_running_ = true;
    should_stop_ = false;

    // Create output directory for visualization artifacts
    std::filesystem::create_directories(config_.export_dir);

    scenario_->setup();

    double total_mlups = 0.0;
    uint64_t measured_steps = 0;
    int export_frame_counter = 0;

    while (!should_stop_ && !scenario_->is_finished() && scenario_->current_step() < config_.max_steps) {
        while (is_paused_ && !should_stop_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        if (should_stop_) break;

        auto t0 = std::chrono::high_resolution_clock::now();

        scenario_->step();
        uint64_t step = scenario_->current_step();
        step_count_ = step;

        double mlups = scenario_->current_mlups();
        if (mlups > 0.0) {
            total_mlups += mlups;
            measured_steps++;
            avg_mlups_ = total_mlups / measured_steps;
        }

        if (step_callback_) {
            step_callback_(step, mlups);
        }

        // Live terminal visualization
        if (config_.live_terminal_render && (step % config_.render_interval == 0)) {
            scenario_->render_terminal(std::cout);
        }

        // Export data frame
        if (step % config_.export_interval == 0) {
            scenario_->export_frame(config_.export_dir, export_frame_counter++);
        }

        // Target FPS throttling
        if (config_.target_fps > 0.0) {
            auto t1 = std::chrono::high_resolution_clock::now();
            double elapsed_sec = std::chrono::duration<double>(t1 - t0).count();
            double frame_budget_sec = 1.0 / config_.target_fps;
            if (elapsed_sec < frame_budget_sec) {
                std::this_thread::sleep_for(
                    std::chrono::duration<double>(frame_budget_sec - elapsed_sec)
                );
            }
        }
    }

    // Final frame and telemetry export
    scenario_->export_frame(config_.export_dir, export_frame_counter);
    scenario_->export_telemetry(config_.export_dir);

    is_running_ = false;
}

void SimulationRunner::start_async() {
    if (is_running_) return;
    worker_thread_ = std::thread([this]() { run(); });
}

void SimulationRunner::pause() {
    is_paused_ = true;
}

void SimulationRunner::resume() {
    is_paused_ = false;
}

void SimulationRunner::stop() {
    should_stop_ = true;
    is_paused_ = false;
}

void SimulationRunner::wait() {
    if (worker_thread_.joinable()) {
        worker_thread_.join();
    }
}

uint64_t SimulationRunner::current_step() const {
    return scenario_ ? scenario_->current_step() : step_count_;
}

} // namespace wake
