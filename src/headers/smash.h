//
// Created by mafu on 1/5/2024.
//

#pragma once
#include "particle.h"
#include <chrono>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>
#define M_PI 3.14159265358979323846
namespace Coal {

    using RapidityRange = std::pair<double, double>;
    using RapidityArray = std::vector<RapidityRange>;

    using ParticleArray      = std::vector<Particle>;
    using MultiParticleArray = std::vector<ParticleArray>;
    using ParticleTypeMap    = std::map<int, ParticleArray>;
    //[pdg,particles]
    using ParticleEventMap = std::map<int, ParticleArray>;
    //[eventID,particles]
    using EventsMap = std::map<int, ParticleTypeMap>;
    //[centrility, EventMap]
    using CentralityMap = std::map<std::pair<int, int>, EventsMap>;

    class PreData {
    public:
        EventsMap AllEvents{};

        PreData() = default;

        static EventsMap readFile(const std::string &filename,
                                  const std::string &mode);

        static void getPtArray(const EventsMap &allEvents, int pdg,
                               const std::string &outputFilename,
                               const RapidityArray &rapidityArray,
                               const std::pair<double, int> &ptBins);

        static CentralityMap getCentralityMap(const EventsMap &allEvents,
                                              const YAML::Node &config);


    private:
        static void readFileSmash(const std::string &filename,
                                  EventsMap &allEvents);

        static void readSmashLine(const std::string &line,
                                  ParticleTypeMap &event);

        static void LoadPDGcode(std::unordered_map<int, int> &pdgChargeMap,
                                const std::string &filename);

        static void
        setParticleCharge(Particle &particle,
                          const std::unordered_map<int, int> &pdgChargeMap);

        static void
        readAmptLine(const std::string &line, ParticleTypeMap &event,
                     const std::unordered_map<int, int> &pdgChargeMap);

        static void readFileAmpt(const std::string &filename,
                                 EventsMap &allEvents);

        static void readLineCluster(const std::string &line,
                                 ParticleTypeMap &events);
        static void readFileCluster(const std::string &filename,
                                 EventsMap &allEvents);

        static void readFileBinary(const std::string &filename,
                                   EventsMap &allEvents);
        static int countChargeParticles(const ParticleEventMap &OneEvent);

        static std::map<int, int>
        calculateMultiplicity(const EventsMap &allEvents);

        static double percentile(const std::vector<int> &multiplicity,
                                 double percent);

        static std::map<std::pair<int, int>, int>
        calculateCentralityBounds(const std::map<int, int> &multiplicity,
                                  const YAML::Node &config);

        static std::map<int, std::pair<int, int>>
        classifyAndCountEvents(const EventsMap &allEvents,
                               const YAML::Node &config);
        static void checkAndCreateOutputDir(const std::string &outputFilename);

        static bool fileExistsInCurrentDir(const std::string &name);
    };

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

}// namespace Coal
