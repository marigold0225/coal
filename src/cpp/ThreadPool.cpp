//
// Created by mafu on 4/4/2024.
//
#include "../headers/ThreadPool.h"


Coal::ThreadPool::ThreadPool(const unsigned int num_threads) : stop_all(false) {
    for (unsigned int i = 0; i < num_threads; ++i) {
        workers.emplace_back(&ThreadPool::workerThread, this);
    }
}

Coal::ThreadPool::~ThreadPool() { stop(); }

void Coal::ThreadPool::enqueueTask(std::function<void()> task) {
    {
        std::lock_guard lock(queue_mutex);
        tasks.push(std::move(task));
        ++active_tasks;
    }
    cv.notify_one();
}

void Coal::ThreadPool::waitAll() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    cv_all_done.wait(lock, [this] { return active_tasks == 0; });
}

void Coal::ThreadPool::stop() {
    {
        std::lock_guard lock(queue_mutex);
        stop_all = true;
    }
    cv.notify_all();

    for (auto &worker: workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

void Coal::ThreadPool::workerThread() {
    while (true) {
        std::function<void()> task;
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            cv.wait(lock, [this] { return !tasks.empty() || stop_all; });

            if (stop_all && tasks.empty()) {
                return;
            }

            task = std::move(tasks.front());
            tasks.pop();
        }
        task();
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            --active_tasks;
            if (active_tasks == 0) {
                cv_all_done.notify_all();
            }
        }
    }
}