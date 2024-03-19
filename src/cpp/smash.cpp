//
// Created by mafu on 1/5/2024.
//
#include "../headers/smash.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <sstream>

void Coal::PreData::readSmashLine(const std::string &line,
                                  ParticleTypeMap &event) {
    std::istringstream stream(line);
    Particle particle{};
    double n_coll = 0.0, form_time = 0.0, xsec_fac = 0.0;
    int id = 0, proc_id_origin = 0, proc_type_origin = 0;
    int pdg_mother1   = 0;
    int pdg_mother2   = 0;
    int baryon_number = 0;

    stream >> particle.t >> particle.x >> particle.y >> particle.z >>
            particle.mass >> particle.p0 >> particle.px >> particle.py >>
            particle.pz >> particle.pdg >> id >> particle.charge >> n_coll >>
            form_time >> xsec_fac >> proc_id_origin >> proc_type_origin >>
            particle.freeze_out_time >> pdg_mother1 >> pdg_mother2 >>
            baryon_number;

    if (!stream.fail()) {
        particle.probability = 1.0;
        particle.getFreezeOutPosition();
        event[particle.pdg].push_back(particle);
    } else {
        std::cerr << "Failed to parse line: " << line << std::endl;
    }
}
void Coal::PreData::LoadPDGcode(std::unordered_map<int, int> &pdgChargeMap,
                                const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: cannot open file " + filename);
    }

    int pdg, charge, spin;
    while (file >> pdg >> charge >> spin) {
        // std::cout << "pdg: " << pdg << " charge: " << charge << std::endl;
        pdgChargeMap[pdg] = charge;
    }

    file.close();
}
void Coal::PreData::setParticleCharge(
        Particle &particle, const std::unordered_map<int, int> &pdgChargeMap) {
    if (const auto it = pdgChargeMap.find(particle.pdg);
        it != pdgChargeMap.end()) {
        particle.charge = it->second;
    } else {
        std::cerr << "PDG code " << particle.pdg << " not found in map."
                  << std::endl;
    }
}

void Coal::PreData::readAmptLine(
        const std::string &line, ParticleTypeMap &event,
        const std::unordered_map<int, int> &pdgChargeMap) {
    std::istringstream stream(line);
    Particle particle{};

    stream >> particle.pdg >> particle.px >> particle.py >> particle.pz >>
            particle.mass >> particle.x >> particle.y >> particle.z >>
            particle.freeze_out_time;
    if (!stream.fail()) {
        particle.probability = 1.0;
        particle.t           = particle.freeze_out_time;
        particle.p0          = std::sqrt(
                particle.px * particle.px + particle.py * particle.py +
                particle.pz * particle.pz + particle.mass * particle.mass);
        setParticleCharge(particle, pdgChargeMap);
        event[particle.pdg].push_back(particle);
    } else {
        std::cerr << "Failed to parse line: " << line << std::endl;
    }
}
void Coal::PreData::readFileAmpt(const std::string &filename,
                                 EventsMap &allEvents) {
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


void Coal::PreData::readFileSmash(const std::string &filename,
                                  EventsMap &allEvents) {
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
                if (line.find("scattering_projectile_target yes") !=
                    std::string::npos) {
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

void Coal::PreData::readLineCluster(const std::string &line,
                                    ParticleTypeMap &events) {
    std::istringstream stream(line);
    Particle particle{};

    stream >> particle.pdg >> particle.px >> particle.py >> particle.pz >>
            particle.p0 >> particle.mass >> particle.x >> particle.y >>
            particle.z >> particle.freeze_out_time >> particle.probability;

    if (!stream.fail()) {
        particle.t = particle.freeze_out_time;
        events[particle.pdg].push_back(particle);
    } else {
        std::cerr << "Failed to parse line: " << line << std::endl;
    }
}
void Coal::PreData::readFileCluster(const std::string &filename,
                                    EventsMap &allEvents) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error: cannot open file " + filename);
    }

    std::string firstLine;
    std::getline(file, firstLine);
    int total_events, total_particles;
    std::cout << firstLine << std::endl;
    sscanf(firstLine.c_str(), "#events:%d size:%d", &total_events,
           &total_particles);
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
        readLineCluster(line, currentEvent);
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


void Coal::PreData::readFileBinary(const std::string &filename,
                                   EventsMap &allEvents) {
    std::ifstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        exit(1);
    }

    int total_events;
    int total_particles;
    file.read(reinterpret_cast<char *>(&total_events), sizeof(total_events));
    file.read(reinterpret_cast<char *>(&total_particles),
              sizeof(total_particles));

    int average_particles = total_particles / total_events;
    int remainder         = total_particles % total_events;

    ParticleTypeMap currentEvent{};
    int eventID       = 0;
    int particleCount = 0;

    while (!file.eof()) {
        Particle particle{};
        file.read(reinterpret_cast<char *>(&particle), sizeof(Particle));
        if (file) {
            currentEvent[particle.pdg].push_back(particle);
            particleCount++;

            if (particleCount >=
                average_particles + (eventID < remainder ? 1 : 0)) {
                allEvents[eventID] = currentEvent;
                currentEvent.clear();
                particleCount = 0;
                eventID++;
            }
        }
    }
    if (!currentEvent.empty()) {
        allEvents[eventID] = currentEvent;
    }
}
int Coal::PreData::countChargeParticles(const ParticleEventMap &OneEvent) {
    int chargeParticleCount = 0;
    for (const auto &particle: OneEvent | std::views::values) {
        for (const auto &p: particle) {
            if (p.charge != 0 && std::abs(p.getRapidity()) < 1.0 &&
                std::abs(p.pdg) > 100) {
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
    if (percent <= 0)
        return multiplicity.front();
    if (percent >= 100)
        return multiplicity.back();

    std::vector<int> sorteDate = multiplicity;
    std::ranges::sort(sorteDate);

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
std::map<std::pair<int, int>, int>
Coal::PreData::calculateCentralityBounds(const std::map<int, int> &multiplicity,
                                         const YAML::Node &config) {
    std::vector<int> multiplicityVector;
    for (const auto &values: multiplicity | std::views::values) {
        multiplicityVector.push_back(values);
    }
    const auto centralityLabels =
            config["General"]["Centrality"]["Ranges"]
                    .as<std::vector<std::pair<int, int>>>();
    std::map<std::pair<int, int>, int> centralityBounds;

    for (const auto &centralityLabel: centralityLabels) {
        centralityBounds[centralityLabel] = static_cast<int>(
                percentile(multiplicityVector, 100 - centralityLabel.second));
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
            if (values >= bound) {
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
                  << " has " << count << " events." << std::endl;
        auto centralityDir = config["Output"]["Path"].as<std::string>();
        centralityDir.append("/")
                .append(std::to_string(label.first))
                .append("-")
                .append(std::to_string(label.second));
        checkAndCreateOutputDir(centralityDir);
    }
    return eventCentrality;
}


Coal::EventsMap Coal::PreData::readFile(const std::string &filename,
                                        const std::string &mode) {
    EventsMap allEvents{};
    if (mode == "cluster") {
        readFileCluster(filename, allEvents);
    } else if (mode == "binary") {
        readFileBinary(filename, allEvents);
    } else if (mode == "smash") {
        readFileSmash(filename, allEvents);
    } else if (mode == "ampt") {
        readFileAmpt(filename, allEvents);
    } else {
        std::cerr << "Error: unknown mode " << mode << std::endl;
        exit(1);
    }
    return allEvents;
}

void Coal::PreData::getPtArray(const EventsMap &allEvents, int pdg,
                               const std::string &outputFilename,
                               const RapidityArray &rapidityArray,
                               const std::pair<double, int> &ptBins) {
    ParticleTypeMap particles;

    for (const auto &particleTypeMap: allEvents | std::views::values) {
        if (particleTypeMap.contains(pdg)) {
            for (const auto &particle: particleTypeMap.at(pdg)) {
                particles[particle.pdg].push_back(particle);
            }
        }
    }
    const int total_events = static_cast<int>(allEvents.size());
    std::cout << "total_events: " << total_events << std::endl;
    std::ofstream output(outputFilename);
    std::map<RapidityRange, std::vector<double>> ptArray;
    std::map<RapidityRange, double> yeildArray;
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(ptBins.second, 0.0);
        yeildArray[rap] = 0.0;
    }

    for (const auto &particle: particles[pdg]) {
        for (const auto &rap: rapidityArray) {
            if (double rapidity = particle.getRapidity();
                rap.first < rapidity && rapidity <= rap.second) {
                const double pt = std::sqrt(particle.px * particle.px +
                                            particle.py * particle.py);
                if (const int index = static_cast<int>(pt / ptBins.first);
                    index < ptBins.second && index >= 0) {
                    ptArray[rap][index] += 1.0 / total_events;
                }
                yeildArray[rap] += 1.0 / total_events;
            }
        }
    }
    for (const auto &rap: rapidityArray) {
        output << "Rapidity range: " << rap.first << "<y<" << rap.second
               << ", cluster yield:"
               << yeildArray[rap] / (rap.second - rap.first) << "\n";
        for (auto i = 0; i < ptBins.second; ++i) {
            const double pt =
                    ptBins.first / 2 + static_cast<double>(i) * ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * ptBins.first *
                                std::abs((rap.second - rap.first)));
            output << pt << " " << ptArray[rap][i] << "\n";
        }
        output << "\n";
    }
    output.close();
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
