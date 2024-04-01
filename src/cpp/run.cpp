//
// Created by mafu on 1/15/2024.
//
#include "../headers/run.h"
#include <iostream>
#include <ranges>

Coal::ConfigParser::ConfigParser(const std::string &filename) {
    YAML::Node config = YAML::LoadFile(filename);
    input             = config;
    general           = config["General"];
    inputFileName     = config["General"]["Path"].as<std::string>() + "/" +
                    config["General"]["Filename"].as<std::string>();
    output        = config["Output"];
    seed          = config["General"]["Seed"].as<int>();
    outputPath    = config["Output"]["Path"].as<std::string>();
    clusterParams = config["ClusterParams"];

    for (const auto &reaction: config["Reactions"]) {
        reactionSwitch[reaction.first.as<std::string>()] =
                reaction.second.as<bool>();
    }
}

void Coal::processReaction(
        const std::string &reactionName, const EventsMap &allEvents,
        const YAML::Node &clusterParamsNode, const ConfigParser &data,
        const std::optional<std::pair<int, int>> &centrality = std::nullopt) {
    const auto cluster = ClusterParams(clusterParamsNode[reactionName]);
    const std::string centralitySuffix =
            centrality ? std::to_string(centrality->first) + "-" +
                                 std::to_string(centrality->second)
                       : "all";

    const auto outputFilename =
            constructFilename(data.outputPath, reactionName, centralitySuffix);
    const auto ptFilename = constructFilename(
            data.outputPath, reactionName + "_pt", centralitySuffix);

    if (data.general["Parallel"]["Enable"].as<bool>()) {
        multithreadedParallelism(allEvents, outputFilename, ptFilename, cluster,
                                 data.input);
    } else {
        singleThreadedParallelism(allEvents, outputFilename, ptFilename,
                                  cluster, data.output);
    }
}


void Coal::handleReactions(const EventsMap &allEvents,
                           const ConfigParser &data) {

    RandomNumber::getInstance().seed(data.general["Seed"].as<int>());

    const auto rapidityRange = data.output["RapidityRange"].as<RapidityArray>();

    auto calculateParticleProperties = [&](const EventsMap &events,
                                           const std::string &particleName,
                                           const int pdg,
                                           const std::string &suffix) {
        const auto outputFilename =
                constructFilename(data.outputPath, particleName, suffix);
        std::cout << "Calculating " << particleName << " for " << suffix
                  << "...\n";
        PreData::getFlow(events, pdg, outputFilename, rapidityRange,
                       {0.2, 10});
    };

    if (data.general["Centrality"]["Enable"].as<bool>()) {
        std::cout << "Centrality calculations being performed..." << std::endl;

        for (auto CentralityMapOpt =
                     PreData::getCentralityMap(allEvents, data.input);
             const auto &[centrality, events]: CentralityMapOpt) {

            auto suffix = std::to_string(centrality.first) + "-" +
                          std::to_string(centrality.second);

            calculateParticleProperties(events, "proton_flow", 2212,
                                        suffix);
        }

    } else {
        calculateParticleProperties(allEvents,"proton_flow_3sub_plus", 2212, "all");
    }

    for (const auto &[reaction, enabled]: data.reactionSwitch) {
        if (!enabled)
            continue;
        std::cout << "Calculating " << reaction << "...\n";
        if (auto CentralityMapOpt = PreData::getCentralityMap(allEvents, data.input);
            data.general["Centrality"]["Enable"].as<bool>()) {
            for (const auto &[centrality, events]: CentralityMapOpt) {
                std::cout << "Processing for centrality " << centrality.first
                          << "-" << centrality.second << std::endl;
                processReaction(reaction, events, data.clusterParams, data,
                                centrality);
            }
        } else {
            processReaction(reaction, allEvents, data.clusterParams, data);
        }
    }
}
