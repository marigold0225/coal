//
// Created by mafu on 4/4/2024.
//

#pragma once
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace Coal {

    class ThreadPool {
    public:
        explicit ThreadPool(unsigned int num_threads);
        ~ThreadPool();

        ThreadPool(const ThreadPool &)            = delete;
        ThreadPool &operator=(const ThreadPool &) = delete;

        ThreadPool(ThreadPool &&)            = delete;
        ThreadPool &operator=(ThreadPool &&) = delete;
        // void enqueueTask(std::function<void(int)> task);
        void enqueueTask(std::function<void()> task);
        void waitAll();
        void stop();

    private:
        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;
        // std::queue<std::pair<std::function<void(int)>, int>> tasks;
        std::mutex queue_mutex;
        std::condition_variable cv;
        std::condition_variable cv_all_done;
        std::size_t active_tasks=0;
        bool stop_all;

        void workerThread();
        // void workerThread(int thread_id);
    };
}// namespace Coal
