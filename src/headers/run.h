//
// Created by mafu on 1/15/2024.
//
#pragma once
#include "config.h"
#include "fileloader.h"
#include "flow.h"
#include <optional>
#include <vector>
namespace Coal {
    //
    // void processReaction(const std::string &reactionName, const EventsMap &allEvents,
    //                      const ConfigParser &config, ResParamsMap &resolution,
    //                      const std::optional<std::pair<int, int>> &centrality);
    void processReaction_v2(const std::string &reactionName, EventsMap &allEvents,

                            const ConfigParser &config, ResParamsMap_vec &resolution,
                            const std::vector<std::string> &FlowName,
                            const std::optional<std::pair<int, int>> &centrality);

    void handleReactions(const EventsMap &allEvents, const ConfigParser &config);
                                          
}// namespace Coal
