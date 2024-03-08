//
// Created by mafu on 1/15/2024.
//
#include "../headers/run.h"
#include <iostream>
#include <ranges>


void Coal::handleCentralityCalculation(const std::string &configFileName) {
    std::cout << "Reading input file..." << std::endl;

    YAML::Node config        = YAML::LoadFile(configFileName);
    auto General             = config["General"];
    auto output              = config["Output"];
    auto ReactionSwitch      = config["Reactions"];
    auto ClusterParams       = config["ClusterParams"];
    const auto inputFilename = General["Path"].as<std::string>() + "/" +
                               General["Filename"].as<std::string>();
    const auto allEvents =
            PreData::readFile(inputFilename, General["Mode"].as<std::string>());
    RandomNumber::getInstance().seed(General["Seed"].as<int>());
    std::cout << "Centrality calculations being performed..." << std::endl;

    const auto CentralityMap = PreData::getCentralityMap(allEvents, config);

    auto rapidityArray = output["RapidityRange"].as<RapidityArray>();

    std::vector<std::string> centralityLabels;
    for (const auto &centrality: CentralityMap | std::views::keys) {
        centralityLabels.emplace_back(std::to_string(centrality.first) + "-" +
                                      std::to_string(centrality.second));
    }
    if (!checkFileExits(output["Path"].as<std::string>(), centralityLabels,
                        "proton_pt")) {
        std::cout << "Calculating proton transverse..." << std::endl;
        for (const auto &centrality: CentralityMap | std::views::keys) {
            const auto &events = CentralityMap.at(centrality);

            std::string outputFilename = constructFilename(
                    output["Path"].as<std::string>(), "proton_pt",
                    std::to_string(centrality.first) + "-" +
                            std::to_string(centrality.second));
            PreData::getPtArray(events, 2212, outputFilename, rapidityArray,
                                {0.2, 10});
        }
    } else {
        std::cout << "proton transverse already calculated..."
                  << "\n";
        std::cout << "start cluster coalesce..."
                  << "\n";
    }
    if (ReactionSwitch["Deuteron"].as<bool>()) {
        prapareEachReaction(config, CentralityMap, "Deuteron");
    }
    if (ReactionSwitch["Helium3"].as<bool>()) {
        prapareEachReaction(config, CentralityMap, "Helium3");
    }
    if (ReactionSwitch["Helium4"].as<bool>()) {
        prapareEachReaction(config, CentralityMap, "Helium4");
    }
    if (ReactionSwitch["Li5"].as<bool>()) {
        prapareEachReaction(config, CentralityMap, "Li5");
    }
    if (ReactionSwitch["Li6"].as<bool>()) {
        prapareEachReaction(config, CentralityMap, "Li6");
    }
    if (ReactionSwitch["Be8"].as<bool>()) {
        prapareEachReaction(config, CentralityMap, "Be8");
    }
}
void Coal::handleNoCentralityCalculation(const std::string &configFileName) {
    YAML::Node config        = YAML::LoadFile(configFileName);
    auto General             = config["General"];
    auto output              = config["Output"];
    auto ReactionSwitch      = config["Reactions"];
    auto ClusterParams       = config["ClusterParams"];
    const auto inputFilename = General["Path"].as<std::string>() + "/" +
                               General["Filename"].as<std::string>();
    RandomNumber::getInstance().seed(General["Seed"].as<int>());
    const auto allEvents =
            PreData::readFile(inputFilename, General["Mode"].as<std::string>());

    std::string protonPTFilename = output["Path"].as<std::string>() + "/" +
                                   "proton_pt" + "all" + ".dat";
    if (!checkFileExits(output["Path"].as<std::string>(), {"all"},
                        "proton_pt")) {
        std::cout << "Calculating proton transverse..."
                  << "\n";
        PreData::getPtArray(allEvents, 2212, protonPTFilename,
                            output["RapidityRange"].as<RapidityArray>(),
                            {0.2, 10});
    } else {
        std::cout << "proton transverse already calculated..."
                  << "\n";
    }

    if (ReactionSwitch["Deuteron"].as<bool>()) {
        std::cout << "Calculating Deuteron..."
                  << "\n";
        std::string outputFilename =
                output["Path"].as<std::string>() + "/" + "Deuteron" + ".dat";
        std::string ptFilename = output["Path"].as<std::string>() + "/" +
                                 "Deuteron" + "_pt" + ".dat";
        Coal::ClusterParams cluster(ClusterParams["Deuteron"]);
        if (General["Parallel"]["Enable"].as<bool>()) {
            clusterThreadV2(allEvents, outputFilename, ptFilename, cluster,
                            config);
        } else {
            clusterMatrix(allEvents, outputFilename, ptFilename, cluster,
                          output);
        }
    }
    if (ReactionSwitch["Helium3"].as<bool>()) {
        std::cout << "Calculating Helium3..."
                  << "\n";
        std::string outputFilename =
                output["Path"].as<std::string>() + "/" + "Helium3" + ".dat";
        std::string ptFilename = output["Path"].as<std::string>() + "/" +
                                 "Helium3" + "_pt" + ".dat";
        Coal::ClusterParams cluster(ClusterParams["Helium3"]);
        if (General["Parallel"]["Enable"].as<bool>()) {
            clusterThreadV2(allEvents, outputFilename, ptFilename, cluster,
                            config);
        } else {
            clusterMatrix(allEvents, outputFilename, ptFilename, cluster,
                          output);
        }
    }
    if (ReactionSwitch["Helium4"].as<bool>()) {
        std::cout << "Calculating Helium4..."
                  << "\n";
        std::string outputFilename =
                output["Path"].as<std::string>() + "/" + "Helium4" + ".dat";
        std::string ptFilename = output["Path"].as<std::string>() + "/" +
                                 "Helium4" + "_pt" + ".dat";
        Coal::ClusterParams cluster(ClusterParams["Helium4"]);
        if (General["Parallel"]["Enable"].as<bool>()) {
            clusterThreadV2(allEvents, outputFilename, ptFilename, cluster,
                            config);
        } else {
            clusterMatrix(allEvents, outputFilename, ptFilename, cluster,
                          output);
        }
    }
    if (ReactionSwitch["Li5"].as<bool>()) {
        std::cout << "Calculating Li5..."
                  << "\n";
        std::string outputFilename =
                output["Path"].as<std::string>() + "/" + "Li5" + ".dat";
        std::string ptFilename =
                output["Path"].as<std::string>() + "/" + "Li5" + "_pt" + ".dat";
        Coal::ClusterParams cluster(ClusterParams["Li5"]);
        if (General["Parallel"]["Enable"].as<bool>()) {
            clusterThreadV2(allEvents, outputFilename, ptFilename, cluster,
                            config);
        } else {
            clusterMatrix(allEvents, outputFilename, ptFilename, cluster,
                          output);
        }
    }
    if (ReactionSwitch["Li6"].as<bool>()) {
        std::cout << "Calculating Li6..."
                  << "\n";
        std::string outputFilename =
                output["Path"].as<std::string>() + "/" + "Li6" + ".dat";
        std::string ptFilename =
                output["Path"].as<std::string>() + "/" + "Li6" + "_pt" + ".dat";
        Coal::ClusterParams cluster(ClusterParams["Li6"]);
        if (General["Parallel"]["Enable"].as<bool>()) {
            clusterThreadV2(allEvents, outputFilename, ptFilename, cluster,
                            config);
        } else {
            clusterMatrix(allEvents, outputFilename, ptFilename, cluster,
                          output);
        }
    }
    if (ReactionSwitch["Be8"].as<bool>()) {
        std::cout << "Calculating Be8..."
                  << "\n";
        std::string outputFilename =
                output["Path"].as<std::string>() + "/" + "Be8" + ".dat";
        std::string ptFilename =
                output["Path"].as<std::string>() + "/" + "Be8" + "_pt" + ".dat";
        Coal::ClusterParams cluster(ClusterParams["Be8"]);
        if (General["Parallel"]["Enable"].as<bool>()) {
            clusterThreadV2(allEvents, outputFilename, ptFilename, cluster,
                            config);
        } else {
            clusterMatrix(allEvents, outputFilename, ptFilename, cluster,
                          output);
        }
    }
}
void Coal::prapareEachReaction(const YAML::Node &config,
                               const CentralityMap &centralityMap,
                               const std::string &reactionName) {
    std::cout << "Calculating " << reactionName << "..."
              << "\n";

    const ClusterParams cluster(config["ClusterParams"][reactionName]);
    for (const auto &centrality: centralityMap | std::views::keys) {
        const auto &events         = centralityMap.at(centrality);
        std::string outputFilename = constructFilename(
                config["Output"]["Path"].as<std::string>(), reactionName + "MC",
                std::to_string(centrality.first) + "-" +
                        std::to_string(centrality.second));
        std::string ptFilename =
                constructFilename(config["Output"]["Path"].as<std::string>(),
                                  reactionName + "MC" + "_pt",
                                  std::to_string(centrality.first) + "-" +
                                          std::to_string(centrality.second));
        if (config["General"]["Parallel"]["Enable"].as<bool>()) {
            clusterThreadV2(events, outputFilename, ptFilename, cluster,
                            config);
        } else {
            clusterMatrix(events, outputFilename, ptFilename, cluster,
                          config["Output"]);
        }
    }
}
void Coal::prapareSpecReaction(const YAML::Node &config,
                               const RapidityArray &rapidityRange,
                               const CentralityMap &centralityMap,
                               const std::string &fromfileName,
                               const std::string &reactionName) {
    CentralityMap newcentralityMap;
    for (const auto &centrality: centralityMap | std::views::keys) {
        auto readFilename = constructFilename(
                config["General"]["Path"].as<std::string>(), fromfileName,
                std::to_string(centrality.first) + "-" +
                        std::to_string(centrality.second));
        newcentralityMap[centrality] = PreData::readFile(readFilename, "spec");
    }
    for (const auto &centrality: newcentralityMap | std::views::keys) {
        const auto &events         = newcentralityMap.at(centrality);
        std::string outputFilename = constructFilename(
                config["General"]["Path"].as<std::string>(), reactionName,
                std::to_string(centrality.first) + "-" +
                        std::to_string(centrality.second));
        std::string ptFilename =
                constructFilename(config["General"]["Path"].as<std::string>(),
                                  reactionName + "_pt",
                                  std::to_string(centrality.first) + "-" +
                                          std::to_string(centrality.second));
        ClusterParams cluster(config["ClusterParams"][reactionName]);
        clusterMix(events, outputFilename, ptFilename, cluster, rapidityRange);
    }
}
