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
        const ConfigParser &data, ResParamsMap &resolution,
        const std::optional<std::pair<int, int>> &centrality = std::nullopt) {

    const auto cluster       = ClusterParams(data.clusterParams[reactionName]);
    const auto rapidityArray = data.rapidityArray;
    const std::string Suffix = centrality ? std::to_string(centrality->first) + "-" +
                                                    std::to_string(centrality->second)
                                          : "0-100";

    const auto clusterFileName = constructFilename(data.outputPath, reactionName, Suffix);

    const auto ptFileName =
            constructFilename(data.outputPath, reactionName + "_flow", Suffix);

    std::optional<unsigned int> max_threads;

    if (data.general["Parallel"]["Enable"].as<bool>()) {
        max_threads = data.general["Parallel"]["Cores"].as<unsigned int>();
    }

    const auto targetParticles =
            ClusterCalculator::Coalescence(allEvents, cluster, resolution, max_threads);

    auto result =
            calculateFlowMatrix(targetParticles, rapidityArray, cluster, resolution);

    const outputInterface outputData("pt");
    outputData.output(result, ptFileName, resolution);
    if (data.output["Extended"].as<bool>()) {
        const outputInterface outputcluster("cluster");
        outputcluster.output(targetParticles, clusterFileName, cluster);
    }
}

void Coal::handleReactions(const EventsMap &allEvents, const ConfigParser &data) {

    RandomNumber::getInstance().seed(data.seed);

    const auto logger = Logger::getInstance().get();

    const auto rapidityArray = data.rapidityArray;

    ResParamsMap resulution = getResolutionMap(allEvents, 1);

    std::optional<CentralityMap> CentralityMapOpt;

    if (data.general["Centrality"]["Enable"].as<bool>()) {
        CentralityMapOpt = Centrality::getCentralityMap(allEvents, data.total);
    }

    CentralityMap processingMap;
    if (CentralityMapOpt) {
        processingMap = *CentralityMapOpt;
    } else {
        processingMap[std::make_pair(0, 100)] = allEvents;
    }

    auto calculateParticleProperties = [&](const EventsMap &events,
                                           const std::string &particleName, const int pdg,
                                           const std::string &suffix) {
        const auto outputFilename =
                constructFilename(data.outputPath, particleName, suffix);
        logger->info("Calculating {} for {}...", particleName, suffix);
        auto result = calculateFlow(events, pdg, rapidityArray, {0.2, 15}, resulution);
        const outputInterface outputData("pt");
        outputData.output(result, outputFilename, resulution);
    };

    for (const auto &[centrality, events]: processingMap) {
        std::string suffix = CentralityMapOpt ? std::to_string(centrality.first) + "-" +
                                                        std::to_string(centrality.second)
                                              : "0-100";
        calculateParticleProperties(events, "proton_flow_long", 2212, suffix);

        for (const auto &[reactionName, enabled]: data.reactionSwitch) {
            if (!enabled)
                continue;
            logger->info("Calculating {} for {}...", reactionName, suffix);
            processReaction(reactionName, events, data, resulution, centrality);
        }
    }
}
