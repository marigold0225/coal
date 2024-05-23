//
// Created by mafu on 4/1/2024.
//
#pragma once
#include "shortcut.h"
#include <map>
#include <string>
#include <unordered_map>

namespace Coal {

    class FileLoader {
    public:
        FileLoader() = default;

        static EventsMap readFile(const std::string &filename, const std::string &mode);

    private:

        static void readMatlabLine(const std::string &line, ParticleTypeMap &event);

        static void readMusicLine(const std::string &line, ParticleTypeMap &event);
        static void readMusicFile(const std::string &filename, EventsMap &allEvents);
        static void readMusicFileCustom(const std::string &filename, EventsMap &allEvents);
        static void readMatlabFile(const std::string &filename, EventsMap &allEvents);

        static void readSmashLine(const std::string &line, ParticleTypeMap &event);

        static void readFileSmash(const std::string &filename, EventsMap &allEvents);

        static void LoadPDGcode(std::unordered_map<int, int> &pdgChargeMap,
                                const std::string &filename);

        static void
        setAmptParticleCharge(Particle &particle,
                              const std::unordered_map<int, int> &pdgChargeMap);

        static void readAmptLine(const std::string &line, ParticleTypeMap &event,
                                 const std::unordered_map<int, int> &pdgChargeMap);

        static void readUrQMDLine(const std::string &line, ParticleTypeMap &event,
                                  const std::unordered_map<int, int> &pdgChargeMap);


        static void readFileUrQMD(const std::string &filename, EventsMap &allEvents);

        static void readFileAmpt(const std::string &filename, EventsMap &allEvents);

        static void readClusterLine(const std::string &line, ParticleTypeMap &event);

        static void readFileCluster(const std::string &filename, EventsMap &allEvents);
    };
}// namespace Coal
