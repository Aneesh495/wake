#pragma once

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

namespace wake {

class ThreadPool {
public:
    explicit ThreadPool(size_t threads = std::thread::hardware_concurrency());
    ~ThreadPool();

    // Enqueue a callable task and return a std::future for result retrieval
    template <class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::invoke_result<F, Args...>::type>;

    // Parallel for loop chunking across worker threads
    void parallel_for(
        int start,
        int end,
        const std::function<void(int start_chunk, int end_chunk)>& chunk_fn
    );

    size_t worker_count() const noexcept { return workers_.size(); }

    void wait_idle();

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;

    std::mutex queue_mutex_;
    std::condition_variable cv_task_;
    std::condition_variable cv_finished_;
    std::atomic<bool> stop_{false};
    std::atomic<int> active_tasks_{0};
};

template <class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
    -> std::future<typename std::invoke_result<F, Args...>::type> {
    using return_type = typename std::invoke_result<F, Args...>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...)
    );

    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        if (stop_) {
            throw std::runtime_error("Cannot enqueue task on stopped ThreadPool");
        }
        tasks_.emplace([task]() { (*task)(); });
    }
    cv_task_.notify_one();
    return res;
}

} // namespace wake
