//
// Created by mafu on 4/2/2024.
//
#pragma once

#include <chrono>
#include <iostream>

namespace Coal {

    using SystemTimePoint = std::chrono::time_point<std::chrono::system_clock>;

    using SystemClock = std::chrono::system_clock;

    using SystemTimeSpan = SystemClock::duration;

    template<typename Duration> void printTime(const Duration &duration) {
        if (auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
            hours.count() > 1) {
            std::cout << "time used: " << hours.count() << " h" << std::endl;
        } else {
            if (auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                minutes.count() > 1) {
                std::cout << "time used: " << minutes.count() << " min" << std::endl;
            } else {
                if (auto seconds =
                            std::chrono::duration_cast<std::chrono::seconds>(duration);
                    seconds.count() > 1) {
                    std::cout << "time used: " << seconds.count() << " s" << std::endl;
                } else {
                    if (auto milliseconds =
                                std::chrono::duration_cast<std::chrono::milliseconds>(
                                        duration);
                        milliseconds.count() > 1) {
                        std::cout << "time used: " << milliseconds.count() << " ms"
                                  << std::endl;
                    }
                }
            }
        }
    }
}// namespace Coal
