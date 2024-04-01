//
// Created by mafu on 4/1/2024.
//
#pragma once
#include "particle.h"
#include <map>
#include <string>
#include <vector>
#include <unordered_map>

namespace Coal {

    using RapidityRange      = std::pair<double, double>;
    using RapidityArray      = std::vector<RapidityRange>;
    using Pair               = std::pair<int, int>;
    using ParticleArray      = std::vector<Particle>;
    using MultiParticleArray = std::vector<ParticleArray>;
    using ParticleTypeMap    = std::map<int, ParticleArray>;
    //[pdg,particles]
    using ParticleEventMap = std::map<int, ParticleArray>;
    //[eventID,particles]
    using EventsMap = std::map<int, ParticleTypeMap>;
    //[centrility, EventMap]
    using CentralityMap = std::map<Pair, EventsMap>;
    //[eventID,EventPlane angle]
    using EventPlaneMap = std::map<int, double>;


    class FileLoader {
    public:
        FileLoader() = default;

        static EventsMap readFile(const std::string &filename, const std::string &mode);

    private:
        static void readSmashLine(const std::string &line,
                                  ParticleTypeMap &event);

        static void readFileSmash(const std::string &filename,
                                  EventsMap &allEvents);

        static void LoadPDGcode(std::unordered_map<int, int> &pdgChargeMap,
                                const std::string &filename);

        static void
        setAmptParticleCharge(Particle &particle,
                          const std::unordered_map<int, int> &pdgChargeMap);

        static void
        readAmptLine(const std::string &line, ParticleTypeMap &event,
                     const std::unordered_map<int, int> &pdgChargeMap);

        static void readFileAmpt(const std::string &filename,
                                 EventsMap &allEvents);

        static void readClusterLine(const std::string &line,
                                    ParticleTypeMap &event);

        static void readFileCluster(const std::string &filename, EventsMap &allEvents);


    };
}

