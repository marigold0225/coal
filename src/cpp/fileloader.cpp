//
// Created by mafu on 4/1/2024.
//
#include "../headers/fileloader.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

Coal::EventsMap Coal::FileLoader::readFile(const std::string &filename,
                                           const std::string &mode) {
    EventsMap allEvents{};
    if (mode == "cluster") {
        readFileCluster(filename, allEvents);
    } else if (mode == "smash") {
        readFileSmash(filename, allEvents);
    } else if (mode == "ampt") {
        readFileAmpt(filename, allEvents);
    } else if (mode == "urqmd") {
        readFileUrQMD(filename, allEvents);
    } else {
        throw std::runtime_error("Error: mode " + mode + " not recognized.");
    }
    return allEvents;
}

void Coal::FileLoader::readSmashLine(const std::string &line, ParticleTypeMap &event) {
    std::istringstream stream(line);
    Particle particle{};
    double n_coll = 0.0, form_time = 0.0, xsec_fac = 0.0;
    int id = 0, proc_id_origin = 0, proc_type_origin = 0;
    int pdg_mother1   = 0;
    int pdg_mother2   = 0;
    int baryon_number = 0;

    stream >> particle.t >> particle.x >> particle.y >> particle.z >> particle.mass >>
            particle.p0 >> particle.px >> particle.py >> particle.pz >> particle.pdg >>
            id >> particle.charge >> n_coll >> form_time >> xsec_fac >> proc_id_origin >>
            proc_type_origin >> particle.freeze_out_time >> pdg_mother1 >> pdg_mother2 >>
            baryon_number;

    if (!stream.fail()) {
        particle.probability = 1.0;
        particle.getFreezeOutPosition();
        event[particle.pdg].push_back(particle);
    } else {
        throw std::runtime_error("Failed to parse line: " + line);
    }
}
void Coal::FileLoader::readFileSmash(const std::string &filename, EventsMap &allEvents) {
    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error("Error: cannot open file " + filename);
    }

    std::string line;
    ParticleTypeMap currentEvent{};
    int currentEventID = 0;
    bool eventStarted  = false;
    bool inEvent       = false;
    bool eventValid    = false;

    while (std::getline(file, line)) {
        if (line.starts_with("# event")) {
            if (line.find("out") != std::string::npos) {
                if (eventStarted && eventValid) {
                    allEvents[currentEventID] = currentEvent;
                }
                currentEvent = ParticleTypeMap{};
                std::istringstream stream(line);
                std::string temp;
                stream >> temp >> temp >> currentEventID;
                eventStarted = true;
                inEvent      = true;
                eventValid   = false;
            } else if (line.find("end") != std::string::npos && inEvent) {
                if (line.find("scattering_projectile_target yes") != std::string::npos) {
                    eventValid = true;
                }
                inEvent = false;
            }
        } else if (inEvent) {
            readSmashLine(line, currentEvent);
        }
    }

    if (eventStarted && eventValid) {
        allEvents[currentEventID] = currentEvent;
    }
}
void Coal::FileLoader::LoadPDGcode(std::unordered_map<int, int> &pdgChargeMap,
                                   const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: cannot open file " + filename);
    }

    int pdg, charge, spin;
    while (file >> pdg >> charge >> spin) {

        pdgChargeMap[pdg] = charge;
    }

    file.close();
}
void Coal::FileLoader::setAmptParticleCharge(
        Particle &particle, const std::unordered_map<int, int> &pdgChargeMap) {

    if (const auto it = pdgChargeMap.find(particle.pdg); it != pdgChargeMap.end()) {
        particle.charge = it->second;
    } else {
        std::cerr << "PDG code " << particle.pdg << " not found in map." << std::endl;
    }
}
void Coal::FileLoader::readAmptLine(const std::string &line, ParticleTypeMap &event,
                                    const std::unordered_map<int, int> &pdgChargeMap) {

    std::istringstream stream(line);
    Particle particle{};

    stream >> particle.pdg >> particle.px >> particle.py >> particle.pz >>
            particle.mass >> particle.x >> particle.y >> particle.z >>
            particle.freeze_out_time;
    if (!stream.fail()) {
        particle.probability = 1.0;
        particle.t           = particle.freeze_out_time;
        particle.p0 =
                std::sqrt(particle.px * particle.px + particle.py * particle.py +
                          particle.pz * particle.pz + particle.mass * particle.mass);
        setAmptParticleCharge(particle, pdgChargeMap);
        event[particle.pdg].push_back(particle);
    } else {
        std::cerr << "Failed to parse line: " << line << std::endl;
    }
}
void Coal::FileLoader::readUrQMDLine(const std::string &line, ParticleTypeMap &event,
                                     const std::unordered_map<int, int> &pdgChargeMap) {

    std::istringstream stream(line);
    Particle particle{};

    int number = 0;
    stream >> number >> particle.pdg >> particle.px >> particle.py >> particle.pz >>
            particle.p0 >> particle.mass >> particle.x >> particle.y >> particle.z >>
            particle.freeze_out_time;
    if (!stream.fail()) {
        particle.probability = 1.0;
        particle.t           = particle.freeze_out_time;
        setAmptParticleCharge(particle, pdgChargeMap);
        event[particle.pdg].push_back(particle);
    } else {
        std::cerr << "Failed to parse line: " << line << std::endl;
    }
}
void Coal::FileLoader::readFileUrQMD(const std::string &filename, EventsMap &allEvents) {

    std::ifstream file(filename);
    std::string line;
    std::unordered_map<int, int> pdgChargeMap;
    LoadPDGcode(pdgChargeMap, "input/AuAuPID.txt");
    ParticleTypeMap currentEvent{};
    int currentEventID      = 0;
    int particleInEvent     = 0;
    int currentPaticleCount = 0;
    double impactParams     = 0.0;
    double other            = 0;

    for (auto i = 0; i < 3; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        if (currentPaticleCount == 0) {
            std::istringstream eventHeader(line);
            eventHeader >> currentEventID >> particleInEvent >> impactParams >> other;
            currentEvent.clear();
        } else {
            readUrQMDLine(line, currentEvent, pdgChargeMap);
        }

        if (++currentPaticleCount > particleInEvent) {
            if (!currentEvent.empty()) {
                allEvents[currentEventID] = currentEvent;
            }
            currentPaticleCount = 0;
        }
    }
    if (!currentEvent.empty()) {
        allEvents[currentEventID] = currentEvent;
    }
}

void Coal::FileLoader::readFileAmpt(const std::string &filename, EventsMap &allEvents) {
    std::ifstream file(filename);
    std::string line;
    std::unordered_map<int, int> pdgChargeMap;
    LoadPDGcode(pdgChargeMap, "input/AuAuPID.txt");
    ParticleTypeMap currentEvent{};
    int currentEventID      = 0;
    int particleInEvent     = 0;
    int currentPaticleCount = 0;
    int other               = 0;

    while (std::getline(file, line)) {
        if (currentPaticleCount == 0) {
            std::istringstream eventHeader(line);
            eventHeader >> currentEventID >> other >> particleInEvent;
            currentEvent.clear();
        } else {
            readAmptLine(line, currentEvent, pdgChargeMap);
        }

        if (++currentPaticleCount > particleInEvent) {
            if (!currentEvent.empty()) {
                allEvents[currentEventID] = currentEvent;
            }
            currentPaticleCount = 0;
        }
    }
    if (!currentEvent.empty()) {
        allEvents[currentEventID] = currentEvent;
    }
}
void Coal::FileLoader::readClusterLine(const std::string &line, ParticleTypeMap &event) {
    std::istringstream stream(line);
    Particle particle{};

    stream >> particle.pdg >> particle.px >> particle.py >> particle.pz >> particle.p0 >>
            particle.mass >> particle.x >> particle.y >> particle.z >>
            particle.freeze_out_time >> particle.probability;

    if (!stream.fail()) {
        particle.t = particle.freeze_out_time;
        event[particle.pdg].push_back(particle);
    } else {
        std::cerr << "Failed to parse line: " << line << std::endl;
    }
}
void Coal::FileLoader::readFileCluster(const std::string &filename,
                                       EventsMap &allEvents) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error: cannot open file " + filename);
    }

    std::string firstLine;
    std::getline(file, firstLine);
    int total_events, total_particles;
    std::cout << firstLine << std::endl;
    sscanf(firstLine.c_str(), "#events:%d size:%d", &total_events, &total_particles);
    int average_particles = total_particles / total_events;
    int remainder         = total_particles % total_events;
    std::string line;
    ParticleTypeMap currentEvent{};
    int eventID   = 0;
    int lineCount = 0;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        readClusterLine(line, currentEvent);
        lineCount++;

        if (lineCount >= average_particles + (eventID < remainder ? 1 : 0)) {
            allEvents[eventID] = currentEvent;
            currentEvent.clear();
            lineCount = 0;
            eventID++;
        }
    }
    if (!currentEvent.empty()) {
        allEvents[eventID] = currentEvent;
    }
}