//
// Created by mafu on 1/15/2024.
//
#pragma once
#include "clusterMix.h"
namespace Coal {

    struct ConfigParser {
        YAML::Node input;
        std::string inputFileName;
        std::string outputPath;
        std::map<std::string, bool> reactionSwitch;
        YAML::Node clusterParams;
        int seed;
        YAML::Node general;
        YAML::Node output;
        explicit ConfigParser(const std::string &filename);
    };

    void processReaction(const std::string &reactionName,
                         const EventsMap &allEvents,
                         const YAML::Node &clusterParamsNode,
                         const ConfigParser &data,
                         const std::optional<std::pair<int, int>> &centrality);

    void handleReactions(const EventsMap &allEvents, const ConfigParser &data);

}// namespace Coal
