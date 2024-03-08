//
// Created by mafu on 1/15/2024.
//
#pragma once
#include "clusterMix.h"
namespace Coal {

    void prapareEachReaction(const YAML::Node &config,
                             const CentralityMap &centralityMap,
                             const std::string &reactionName);
    void prapareSpecReaction(const YAML::Node &config,
                             const RapidityArray &rapidityRange,
                             const CentralityMap &centralityMap,
                             const std::string &fromfileName,
                             const std::string &reactionName);

    void handleCentralityCalculation(const std::string &configFileName);

    void handleNoCentralityCalculation(const std::string &configFileName);
}// namespace Coal
