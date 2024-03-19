//
// Created by mafu on 1/15/2024.
//
#pragma once
#include "clusterMix.h"
namespace Coal {

    struct ConfigData {
        YAML::Node input;
        std::string inputFileName;
        std::string outputPath;
        std::map<std::string, bool> reactionSwitch;
        YAML::Node clusterParams;
        int seed;
        YAML::Node general;
        YAML::Node output;
    };

    ConfigData readCondig(const std::string &fileName);

    void processReaction(const std::string &reactionName,
                         const EventsMap &allEvents,
                         const YAML::Node &clusterParamsNode,
                         const ConfigData &data,
                         const std::optional<std::pair<int, int>> &centrality);

    void handleReactions(const EventsMap &allEvents, const ConfigData &data);

    void prapareEachReaction(const YAML::Node &config,
                             const CentralityMap &centralityMap,
                             const std::string &reactionName);
    void prapareSpecReaction(const YAML::Node &config,
                             const CentralityMap &centralityMap,
                             const std::string &fromfileName,
                             const std::string &reactionName);

    void handleCentralityCalculation(const std::string &configFileName);

    void handleNoCentralityCalculation(const std::string &configFileName);
}// namespace Coal
