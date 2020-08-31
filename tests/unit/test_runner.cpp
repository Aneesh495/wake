#include "wake/runner/thread_pool.hpp"
#include "wake/runner/simulation_runner.hpp"
#include "wake/scenarios/karman_vortex.hpp"
#include <iostream>
#include <cassert>
#include <atomic>

using namespace wake;

void test_thread_pool() {
    std::cout << "Testing ThreadPool task enqueue and parallel_for..." << std::endl;

    ThreadPool pool(4);
    assert(pool.worker_count() == 4);

    // Enqueue future
    auto fut = pool.enqueue([](int a, int b) { return a + b; }, 19, 23);
    assert(fut.get() == 42);

    // Parallel for
    const int N = 1000;
    std::vector<int> data(N, 0);

    pool.parallel_for(0, N, [&data](int start, int end) {
        for (int i = start; i < end; ++i) {
            data[i] = i * 2;
        }
    });

    for (int i = 0; i < N; ++i) {
        assert(data[i] == i * 2);
    }

    std::cout << "ThreadPool PASSED." << std::endl;
}

void test_simulation_runner() {
    std::cout << "Testing SimulationRunner with Karman scenario..." << std::endl;

    auto scenario = std::make_shared<KarmanVortexScenario>(60, 30, 100.0f, 0.05f, 8.0f);

    RunnerConfig config;
    config.max_steps = 30;
    config.render_interval = 1000; // disable during unit test
    config.export_interval = 1000;
    config.live_terminal_render = false;
    config.export_dir = "./test_run_output";
    config.target_fps = 0.0; // unthrottled

    SimulationRunner runner(scenario, config);

    std::atomic<int> callbacks_received{0};
    runner.set_step_callback([&callbacks_received](uint64_t step, double mlups) {
        (void)step; (void)mlups;
        callbacks_received++;
    });

    runner.run();

    assert(runner.current_step() == 30);
    assert(callbacks_received == 30);
    assert(runner.average_mlups() > 0.0);

    std::filesystem::remove_all("./test_run_output");

    std::cout << "SimulationRunner PASSED (average MLUPS = " << runner.average_mlups() << ")." << std::endl;
}

int main() {
    test_thread_pool();
    test_simulation_runner();
    std::cout << "All runner tests passed successfully!" << std::endl;
    return 0;
}
