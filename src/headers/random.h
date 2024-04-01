//
// Created by mafu on 4/1/2024.
//
#pragma once
#include <random>
#include <chrono>
namespace Coal {
    class RandomNumber {
    public:
        // void seed(const unsigned int seed) { gen.seed(seed); }
        void seed(const int seed) {
            if (seed == -1) {
                const auto now = std::chrono::system_clock::now();
                const auto milliseconds =
                        std::chrono::duration_cast<std::chrono::milliseconds>(
                                now.time_since_epoch())
                                .count();
                gen.seed(static_cast<unsigned int>(milliseconds % 100000000));
            } else {
                gen.seed(static_cast<unsigned int>(seed));
            }
        }

        double uniform(const double min, const double max) {
            std::uniform_real_distribution dis(min, max);
            return dis(gen);
        }

        static RandomNumber &getInstance() {
            static RandomNumber instance;
            return instance;
        }

        RandomNumber(RandomNumber const &) = delete;

        void operator=(RandomNumber const &) = delete;

        std::mt19937 &getGenerator() { return gen; }

    private:
        std::mt19937 gen;

        RandomNumber() {
            std::random_device rd;
            gen = std::mt19937(rd());
        }
    };
}

