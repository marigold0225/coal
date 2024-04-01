//
// Created by mafu on 1/5/2024.
//
#include "../headers/smash.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <oneapi/tbb/task.h>
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
                const double pt = particle.getPT();
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
Coal::ResParams Coal::PreData::getResolution(const ParticleEventMap &OneEvent,
                                             const int N) {

    ResParams v1Params{};
    double Q_x   = 0.0;
    double Q_y   = 0.0;
    double Q_x_A = 0.0;
    double Q_y_A = 0.0;
    double Q_x_B = 0.0;
    double Q_y_B = 0.0;
    double Q_x_C = 0.0;
    double Q_y_C = 0.0;

    for (const auto &particles: OneEvent | std::views::values) {
        for (const auto &particle: particles) {
            const double pseudoRapidity = particle.getPseudoRapidity();
            const double pT             = particle.getPT();
            const double Phi            = particle.getPhi();

            if (pseudoRapidity > 2.81 && pseudoRapidity < 5.09) {
                Q_x += pT * std::cos(Phi);
                Q_y += pT * std::sin(Phi);
            }
            if (pseudoRapidity > 2.81 && pseudoRapidity < 5.09) {
                Q_x_A += pT * std::cos(Phi);
                Q_y_A += pT * std::sin(Phi);
            }
            if (pseudoRapidity > 2.14 && pseudoRapidity < 2.81) {
                Q_x_B += pT * std::cos(Phi);
                Q_y_B += pT * std::sin(Phi);
            }
            if (pseudoRapidity > 0.1 && pseudoRapidity < 2.0) {
                Q_x_C += pT * std::cos(Phi);
                Q_y_C += pT * std::sin(Phi);
            }
        }
    }
    v1Params.Q_x         = Q_x;
    v1Params.Q_y         = Q_y;
    v1Params.Psi_1       = Psi(Q_y, Q_x);
    v1Params.Q_x_A       = Q_x_A;
    v1Params.Q_y_A       = Q_y_A;
    v1Params.Q_x_B       = Q_x_B;
    v1Params.Q_y_B       = Q_y_B;
    const double Phi_1   = Psi(Q_y_A, Q_x_A);
    const double Phi_CD  = Psi(Q_y_B, Q_x_B);
    const double Phi_TPC = Psi(Q_y_C, Q_x_C);

    v1Params.cos_N_AB = std::cos(N * (Phi_1 - Phi_CD));
    v1Params.cos_N_AC = std::cos(N * (Phi_1 - Phi_TPC));
    v1Params.cos_N_BC = std::cos(N * (Phi_CD - Phi_TPC));
    return v1Params;
}

void Coal::PreData::rotateEventPlane(ParticleTypeMap &OneEvent,
                                     const double angle) {
    // Rotate the event plane by the specified angle
    for (auto &particles: OneEvent | std::views::values) {
        for (auto &particle: particles) {
            const double px = particle.px;
            const double py = particle.py;
            const double x  = particle.x;
            const double y  = particle.y;
            particle.px     = px * std::cos(angle) - py * std::sin(angle);
            particle.py     = px * std::sin(angle) + py * std::cos(angle);
            particle.x      = x * std::cos(angle) - y * std::sin(angle);
            particle.y      = x * std::sin(angle) + y * std::cos(angle);
        }
    }
}
void Coal::PreData::getFlow(const EventsMap &allEvents, int pdg,
                          const std::string &outputFilename,
                          const RapidityArray &rapidityArray,
                          const std::pair<double, int> &ptBins) {

    double cos_N_AB = 0.0;
    double cos_N_AC = 0.0;
    double cos_N_BC = 0.0;
    const int total_events = static_cast<int>(allEvents.size());
    std::pair ptRange      = {0.4, 2.0};
    std::map<int, ResParams> ResParamsMap;
    std::map<RapidityRange, std::vector<double>> ptArray;
    std::map<RapidityRange, double> yeildArray;
    std::map<RapidityRange, std::vector<double>> v1Array;
    std::map<RapidityRange, std::vector<double>> v2Array;
    std::map<RapidityRange, std::vector<int>> v1Count;
    std::map<RapidityRange, std::vector<int>> v2Count;
    std::map<RapidityRange, double> v1overpt;
    std::map<RapidityRange, int> v1overptCount;
    std::map<RapidityRange, double> v2overpt;
    std::map<RapidityRange, int> v2overptCount;
    for (const auto &rap: rapidityArray) {
        ptArray[rap]       = std::vector(ptBins.second, 0.0);
        yeildArray[rap]    = 0.0;
        v1Array[rap]       = std::vector(ptBins.second, 0.0);
        v1Count[rap]       = std::vector(ptBins.second, 0);
        v2Array[rap]       = std::vector(ptBins.second, 0.0);
        v2Count[rap]       = std::vector(ptBins.second, 0);
        v1overpt[rap]      = 0.0;
        v1overptCount[rap] = 0;
        v2overpt[rap]      = 0.0;
        v2overptCount[rap] = 0;
    }

    std::ofstream output(outputFilename);

    for (const auto &[event_ID, oneEvent]: allEvents) {
        constexpr int N        = 1;
        ResParamsMap[event_ID] = getResolution(oneEvent, N);
        cos_N_AB += ResParamsMap[event_ID].cos_N_AB;
        cos_N_AC += ResParamsMap[event_ID].cos_N_AC;
        cos_N_BC += ResParamsMap[event_ID].cos_N_BC;
    }
    double cos_N_AB_mean = cos_N_AB / total_events;
    double cos_N_AC_mean = cos_N_AC / total_events;
    double cos_N_BC_mean = cos_N_BC / total_events;

    for (const auto &[eventId, OneEvent]: allEvents) {
        for (const auto &[particleType, particles]: OneEvent) {
            if (particleType == pdg) {
                for (const auto &particle: particles) {
                    for (const auto &rap: rapidityArray) {
                        if (double rapidity = particle.getRapidity();
                            rap.first < rapidity && rapidity <= rap.second) {
                            const double pT  = particle.getPT();
                            const double phi = particle.getPhi();
                            const double psi    = ResParamsMap[eventId].Psi_1;
                            const double v1_obs = std::cos(phi - psi);
                            const double v2_obs = std::cos(2 * (phi - psi));

                            if (const int index =
                                        static_cast<int>(pT / ptBins.first);
                                index < ptBins.second && index >= 0) {
                                ptArray[rap][index] += 1.0 / total_events;
                                v1Array[rap][index] += v1_obs;
                                v1Count[rap][index]++;
                                v2Array[rap][index] += v2_obs;
                                v2Count[rap][index]++;
                            }
                            if (pT >= ptRange.first && pT <= ptRange.second) {
                                v1overpt[rap] += v1_obs;
                                v1overptCount[rap]++;
                                v2overpt[rap] += v2_obs;
                                v2overptCount[rap]++;
                            }
                            yeildArray[rap] += 1.0 / total_events;
                        }
                    }
                }
            }
        }
    }

    const double Ref_v1 =
            std::sqrt(cos_N_AB_mean * cos_N_AC_mean / cos_N_BC_mean);
    const double x      = solveEquation(Ref_v1, 1);
    const double Ref_v2 = getRef(x, 2);

    std::cout << "Ref_v1:" << Ref_v1 << "Ref_v2:" << Ref_v2 << std::endl;

    for (const auto &rap: rapidityArray) {
        output << "Rapidity range: " << rap.first << "<y<" << rap.second
               << ", cluster yield:"
               << yeildArray[rap] / (rap.second - rap.first) << "\n";
        output << "flow over pt range: " << ptRange.first << "<pt<"
               << ptRange.second
               << ", v1:" << v1overpt[rap] / v1overptCount[rap] / Ref_v1
               << ", v2:" << v2overpt[rap] / v2overptCount[rap] / Ref_v2
               << "\n";
        for (auto i = 0; i < ptBins.second; ++i) {
            const double pt =
                    ptBins.first / 2 + static_cast<double>(i) * ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * ptBins.first *
                                std::abs((rap.second - rap.first)));

            v1Array[rap][i] =
                    v1Count[rap][i] > 0
                            ? v1Array[rap][i] / v1Count[rap][i] / Ref_v1
                            : 0.0;

            v2Array[rap][i] =
                    v2Count[rap][i] > 0
                            ? v2Array[rap][i] / v2Count[rap][i] / Ref_v2
                            : 0.0;
            output << pt << " " << ptArray[rap][i] << " " << v1Array[rap][i]
                   << " " << v2Array[rap][i] << "\n";
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
double Coal::PreData::Psi(const double y, const double x) {
    double angle = atan2(y, x);
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    return angle;
}
