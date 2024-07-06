//
// Created by mafu on 1/5/2024.
//

#pragma once
#include "shortcut.h"
#include <map>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>
namespace Coal {

    class Centrality {
    public:
        Centrality() = default;

        static CentralityMap getCentralityMap(const EventsMap &allEvents,
                                              const YAML::Node &config);

    private:

        static int countChargeParticles(const ParticleEventMap &OneEvent);

        static std::map<int, int>
        calculateMultiplicity(const EventsMap &allEvents);

        static double percentile(const std::vector<int> &multiplicity,
                                 double percent);

        static std::map<Pair, Pair>
        calculateCentralityBounds(const std::map<int, int> &multiplicity,
                                  const YAML::Node &config);

        static std::map<int, std::pair<int, int>>
        classifyAndCountEvents(const EventsMap &allEvents,
                               const YAML::Node &config);
        static void checkAndCreateOutputDir(const std::string &outputFilename);

        static bool fileExistsInCurrentDir(const std::string &name);
    };
}// namespace Coal
