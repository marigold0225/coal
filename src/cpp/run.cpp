//
// Created by mafu on 1/15/2024.
//
#include "../headers/run.h"
#include "../headers/Centrality.h"
#include "../headers/ClusterCalculator.h"
#include "../headers/logger.h"
#include "../headers/output.h"
#include "../headers/random.h"

void Coal::processReaction(
        const std::string &reactionName, const EventsMap &allEvents,
        const ConfigParser &config, ResParamsMap &resolution,
        const std::optional<std::pair<int, int>> &centrality = std::nullopt) {

    const auto clusterConfig = ClusterParams(config.clusterParams[reactionName]);
    const auto rapidityArray = config.rapidityArray;
    const std::string Suffix = centrality ? std::to_string(centrality->first) + "-" +
                                                    std::to_string(centrality->second)
                                          : "0-100";

    const auto clusterDir = constructFilename(config.outputPath, reactionName, Suffix);

    const auto ptDir = constructFilename(config.outputPath, reactionName + "_flow", Suffix);

    std::optional<unsigned int> max_threads;

    if (config.general["Parallel"]["Enable"].as<bool>()) {
        max_threads = config.general["Parallel"]["Cores"].as<unsigned int>();
    }

    const auto clusterResult = ClusterCalculator::Coalescence(allEvents, clusterConfig,
                                                              resolution, max_threads);

    auto result =
            calculateFlowMatrix(clusterResult, rapidityArray, clusterConfig, resolution);

    const outputInterface outputData("pt");
    outputData.output(result, ptDir, resolution);
    if (config.output["Extended"].as<bool>()) {
        const outputInterface outputcluster("cluster");
        outputcluster.output(clusterResult, clusterDir, clusterConfig);
    }
}

void Coal::handleReactions(const EventsMap &allEvents, const ConfigParser &config) {

    RandomNumber::getInstance().seed(config.seed);

    const auto logger = Logger::getInstance().get();

    const auto rapidityArray = config.rapidityArray;

    ResParamsMap resulution = getResolutionMap(allEvents, 2);

    std::optional<CentralityMap> CentralityMapOpt;

    if (config.general["Centrality"]["Enable"].as<bool>()) {
        CentralityMapOpt = Centrality::getCentralityMap(allEvents, config.total);
    }

    CentralityMap processingMap;
    if (CentralityMapOpt) {
        processingMap = *CentralityMapOpt;
    } else {
        processingMap[std::make_pair(0, 100)] = allEvents;
    }

    auto calculateParticleProperties = [&](const EventsMap &events,
                                           const std::string &flowName,
                                           const std::string &particleName, const int pdg,
                                           const std::string &suffix) {
        const auto outputDir = constructFilename(config.outputPath, particleName, suffix);
        const auto outputFlowDir = constructFilename(config.outputPath, flowName, suffix);
        logger->info("Calculating {} for {}...", particleName, suffix);
        auto result = calculateFlow(events, pdg, rapidityArray, {0.3, 20}, resulution);
        const outputInterface outputData("pt");
        outputData.output(result, outputFlowDir, resulution);
        // const outputInterface outputParticle("particle");
        // outputParticle.output(events, outputDir, pdg);
    };

    for (const auto &[centrality, events]: processingMap) {
        std::string suffix = CentralityMapOpt ? std::to_string(centrality.first) + "-" +
                                                        std::to_string(centrality.second)
                                              : "0-100";
        calculateParticleProperties(events, "proton_flow", "proton", 2212, suffix);
        // calculateParticleProperties(events, "anti_proton_flow", "antiproton", -2212, suffix);
        // calculateParticleProperties(events, "K_flow", "K", 321, suffix);
        // calculateParticleProperties(events, "anti_K_flow", "antiK", -321, suffix);
        // calculateParticleProperties(events, "Pi_flow", "Pi", 211, suffix);
        // calculateParticleProperties(events, "anti_Pi_flow", "antiPi", -211, suffix);


        for (const auto &[reactionName, enabled]: config.reactionSwitch) {
            if (!enabled)
                continue;
            logger->info("Calculating {} for {}...", reactionName, suffix);
            processReaction(reactionName, events, config, resulution, centrality);
        }
    }
}
