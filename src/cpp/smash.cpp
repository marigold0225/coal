//
// Created by mafu on 1/5/2024.
//
#include "../headers/smash.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
int Coal::PreData::countChargeParticles(const ParticleEventMap &OneEvent) {
    int chargeParticleCount = 0;
    for (const auto &particle: OneEvent | std::views::values) {
        for (const auto &p: particle) {
            if (p.charge != 0 && p.getPseudoRapidity() < 0.0 &&
                p.getPseudoRapidity() > -2.0 && std::abs(p.pdg) > 100) {
                chargeParticleCount++;
            }
        }
    }
    return chargeParticleCount;
}
std::map<int, int>
Coal::PreData::calculateMultiplicity(const EventsMap &allEvents) {
    std::map<int, int> multiplicity;
    for (const auto &[eventID, event]: allEvents) {
        multiplicity[eventID] = countChargeParticles(event);
    }
    return multiplicity;
}
double Coal::PreData::percentile(const std::vector<int> &multiplicity,
                                 const double percent) {
    if (multiplicity.empty()) {
        return 0.0;
    }
    std::vector<int> sorteDate = multiplicity;
    std::ranges::sort(sorteDate);

    if (percent <= 0)
        return sorteDate.front();
    if (percent >= 100)
        return sorteDate.back();


    const double idx =
            static_cast<double>(multiplicity.size() - 1) * percent / 100.0;
    const int idxlower = static_cast<int>(idx);
    const int idxhigh  = idxlower + 1;
    const double frac  = idx - idxlower;

    if (idxhigh >= multiplicity.size()) {
        return sorteDate[idxlower];
    }
    return sorteDate[idxlower] * (1.0 - frac) + sorteDate[idxhigh] * frac;
}
std::map<Coal::Pair, Coal::Pair>
Coal::PreData::calculateCentralityBounds(const std::map<int, int> &multiplicity,
                                         const YAML::Node &config) {
    std::vector<int> multiplicityVector;
    for (const auto &values: multiplicity | std::views::values) {
        multiplicityVector.push_back(values);
    }
    const auto centralityLabels =
            config["General"]["Centrality"]["Ranges"]
                    .as<std::vector<std::pair<int, int>>>();
    std::map<Pair, Pair> centralityBounds;


    for (const auto &centralityLabel: centralityLabels) {

        auto lowerBounds = static_cast<int>(
                percentile(multiplicityVector, 100 - centralityLabel.first));
        auto upperBounds = static_cast<int>(
                percentile(multiplicityVector, 100 - centralityLabel.second));

        centralityBounds[centralityLabel] = {lowerBounds, upperBounds};
    }
    return centralityBounds;
}
void Coal::PreData::checkAndCreateOutputDir(const std::string &outputFilename) {
    if (!std::filesystem::exists(outputFilename)) {
        std::filesystem::create_directories(outputFilename);
    }
}
bool Coal::PreData::fileExistsInCurrentDir(const std::string &name) {
    return exists(std::filesystem::current_path() / name);
}
std::map<int, std::pair<int, int>>
Coal::PreData::classifyAndCountEvents(const EventsMap &allEvents,
                                      const YAML::Node &config) {
    auto multiplicity     = calculateMultiplicity(allEvents);
    auto centralityBounds = calculateCentralityBounds(multiplicity, config);
    std::map<int, std::pair<int, int>> eventCentrality;
    std::map<std::pair<int, int>, int> centralityEventCounts;
    for (const auto &[eventID, values]: multiplicity) {
        bool classified = false;
        for (const auto &[label, bound]: centralityBounds) {
            if (values >= bound.second && values <= bound.first) {
                eventCentrality[eventID] = label;
                centralityEventCounts[label]++;
                classified = true;
                break;
            }
        }
        if (!classified) {
            eventCentrality[eventID] = {-1, -1};
            centralityEventCounts[{-1, -1}]++;
        }
    }
    for (const auto &[label, count]: centralityEventCounts) {
        std::cout << "Centrality " << label.first << "-" << label.second
                  << " has " << count << "events." << std::endl;
        auto centralityDir = config["Output"]["Path"].as<std::string>();
        centralityDir.append("/")
                .append(std::to_string(label.first))
                .append("-")
                .append(std::to_string(label.second));
        checkAndCreateOutputDir(centralityDir);
    }
    return eventCentrality;
}

Coal::CentralityMap Coal::PreData::getCentralityMap(const EventsMap &allEvents,
                                                    const YAML::Node &config) {
    CentralityMap centralityMap;
    auto eventCentrality = classifyAndCountEvents(allEvents, config);
    for (const auto &[eventId, OneEvent]: allEvents) {
        const auto &centralityLabel = eventCentrality[eventId];

        if (centralityLabel.first == -1) {
            continue;
        }
        centralityMap[centralityLabel][eventId] = OneEvent;
    }
    return centralityMap;
}

