//
// Created by mafu on 24-3-26.
//
#include "../headers/flow.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <ranges>
#include <vector>

Coal::flowResult::flowResult(const RapidityArray &rapidityArray, const int bins)
{
    for (const auto &rap: rapidityArray) {
        ptArray[rap]       = std::vector(bins, 0.0);
        yeildArray[rap]    = 0.0;
        v1Array[rap]       = std::vector(bins, 0.0);
        v1Count[rap]       = std::vector(bins, 0);
        v2Array[rap]       = std::vector(bins, 0.0);
        v2Count[rap]       = std::vector(bins, 0);
        v1overpt[rap]      = 0.0;
        v1overptCount[rap] = 0;
        v2overpt[rap]      = 0.0;
        v2overptCount[rap] = 0;
    }
}

Coal::flowResult_vec::flowResult_vec(const RapidityArray &rapidityArray, const int bins,
                                     const std::vector<std::string> &FlowName)
{
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(bins, 0.0);
        yeildArray[rap] = 0.0;
        for (const auto &name: FlowName) {
            flowArray[name][rap]       = std::vector(bins, 0.0);
            flowCount[name][rap]       = std::vector(bins, 0);
            flowOverPt[name][rap]      = 0.0;
            flowOverPtCount[name][rap] = 0;
        }
    }
}

auto Coal::eventPlaneEquation_std(const double x_val, const params *p, double &y) -> int
{
    const double a   = p->a;
    const int    k   = p->k;
    const double arg = x_val * x_val / 4;

    const double I_k_minus_1_over_2 = std::cyl_bessel_i((k - 1) / 2.0, arg);
    const double I_k_plus_1_over_2  = std::cyl_bessel_i((k + 1) / 2.0, arg);

    y = M_PI / 2 / sqrt(2) * x_val * exp(-arg) * (I_k_minus_1_over_2 + I_k_plus_1_over_2) - a;
    return 0;
}

auto Coal::solveEquation_std(const double a, const int k) -> double
{
    const params  par      = {a, k};
    double        x        = 0.5;
    constexpr int max_iter = 1000;
    for (int iter = 0; iter < max_iter; ++iter) {
        double y;
        eventPlaneEquation_std(x, &par, y);
        if (constexpr double tol = 1e-7; std::abs(y) < tol) {
            break;
        }
        constexpr double h = 1e-5;
        double           y1;
        double           y2;
        eventPlaneEquation_std(x + h, &par, y1);
        eventPlaneEquation_std(x - h, &par, y2);
        const double dy = (y1 - y2) / (2 * h);
        x -= y / dy;
    }
    return x;
}

auto Coal::getRef_std(const double x, const int k) -> double
{
    const double arg                = x * x / 4;
    const double I_k_minus_1_over_2 = std::cyl_bessel_i((k - 1) / 2.0, arg);
    const double I_k_plus_1_over_2  = std::cyl_bessel_i((k + 1) / 2.0, arg);
    return M_PI / 2 / sqrt(2) * x * std::exp(-arg) * (I_k_minus_1_over_2 + I_k_plus_1_over_2);
}

auto Coal::Psi(const double y, const double x) -> double
{
    double angle = atan2(y, x);
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    return angle;
}

auto Coal::getResolution(const ParticleEventMap &OneEvent, const int N) -> Coal::ResParams
{

    ResParams params{};
    double    Q_x   = 0.0;
    double    Q_y   = 0.0;
    double    Q_x_A = 0.0;
    double    Q_y_A = 0.0;
    double    Q_x_B = 0.0;
    double    Q_y_B = 0.0;
    double    Q_x_C = 0.0;
    double    Q_y_C = 0.0;

    for (const auto &particles: OneEvent | std::views::values) {
        for (const auto &particle: particles) {
            const double pseudoRapidity = particle.getPseudoRapidity();
            // const double rapidity       = particle.getRapidity();
            const double pT = particle.getPT();
            // const double Phi = particle.getPhi();
            const double Phi = Psi(particle.py, particle.px);
            //AuAu--3GeV--star
            // if (pseudoRapidity > 2.81 && pseudoRapidity < 5.39) {
            //     Q_x += pT * std::cos(Phi);
            //     Q_y += pT * std::sin(Phi);
            // }
            // if (pseudoRapidity > 2.81 && pseudoRapidity < 5.09) {
            //     Q_x_A += pT * std::cos(Phi);
            //     Q_y_A += pT * std::sin(Phi);
            // }
            // if (pseudoRapidity > 2.14 && pseudoRapidity < 2.81) {
            //     Q_x_B += pT * std::cos(Phi);
            //     Q_y_B += pT * std::sin(Phi);
            // }
            // if (pseudoRapidity > 0.1 && pseudoRapidity < 2.0) {
            //     Q_x_C += pT * std::cos(Phi);
            //     Q_y_C += pT * std::sin(Phi);
            // }
            //PbPb--5.02TeV--Alice
            // if (pseudoRapidity > -0.8 && pseudoRapidity < 0.8 && pT > 0.2 && pT < 20) {
            if (pseudoRapidity > -3.7 && pseudoRapidity < -1.7 ||
                pseudoRapidity > 2.8 && pseudoRapidity < 5.1) {
                Q_x += pT * std::cos(N * Phi);
                Q_y += pT * std::sin(N * Phi);
            }
            if (pseudoRapidity > -3.7 && pseudoRapidity < -1.7 ||
                pseudoRapidity > 2.8 && pseudoRapidity < 5.1) {
                Q_x_A += pT * std::cos(N * Phi);
                Q_y_A += pT * std::sin(N * Phi);
            }
            if (pseudoRapidity < 0.8 && pseudoRapidity > 0.0 && pT > 0.2 && pT < 20) {
                Q_x_B += pT * std::cos(N * Phi);
                Q_y_B += pT * std::sin(N * Phi);
            }
            if (pseudoRapidity > -0.8 && pseudoRapidity < 0.0 && pT > 0.2 && pT < 20) {
                Q_x_C += pT * std::cos(N * Phi);
                Q_y_C += pT * std::sin(N * Phi);
            }
        }
    }
    params.Q_x           = Q_x;
    params.Q_y           = Q_y;
    params.Psi_1         = Psi(Q_y, Q_x) / N;
    params.Q_x_A         = Q_x_A;
    params.Q_y_A         = Q_y_A;
    params.Q_x_B         = Q_x_B;
    params.Q_y_B         = Q_y_B;
    params.Q_x_C         = Q_x_C;
    params.Q_y_C         = Q_y_C;
    const double Phi_1   = Psi(Q_y_A, Q_x_A) / N;
    const double Phi_CD  = Psi(Q_y_B, Q_x_B) / N;
    const double Phi_TPC = Psi(Q_y_C, Q_x_C) / N;

    params.cos_N_AB = std::cos(N * (Phi_1 - Phi_CD));
    params.cos_N_AC = std::cos(N * (Phi_1 - Phi_TPC));
    params.cos_N_BC = std::cos(N * (Phi_CD - Phi_TPC));
    return params;
}

auto Coal::getResolution_vec(const ParticleEventMap         &OneEvent,
                             const std::vector<std::string> &N_list) -> Coal::ResParams_vec
{

    ResParams_vec params_vec{};

    for (const auto &N: N_list) {
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
                const double Phi            = Psi(particle.py, particle.px);

                if (pseudoRapidity > -3.7 && pseudoRapidity < -1.7 ||
                    pseudoRapidity > 2.8 && pseudoRapidity < 5.1) {
                    Q_x += pT * std::cos(mapFLowToInt[N] * Phi);
                    Q_y += pT * std::sin(mapFLowToInt[N] * Phi);
                }
                if (pseudoRapidity > -3.7 && pseudoRapidity < -1.7 ||
                    pseudoRapidity > 2.8 && pseudoRapidity < 5.1) {
                    Q_x_A += pT * std::cos(mapFLowToInt[N] * Phi);
                    Q_y_A += pT * std::sin(mapFLowToInt[N] * Phi);
                }
                if (pseudoRapidity < 0.8 && pseudoRapidity > 0.0 && pT > 0.2 && pT < 20) {
                    Q_x_B += pT * std::cos(mapFLowToInt[N] * Phi);
                    Q_y_B += pT * std::sin(mapFLowToInt[N] * Phi);
                }
                if (pseudoRapidity > -0.8 && pseudoRapidity < 0.0 && pT > 0.2 && pT < 20) {
                    Q_x_C += pT * std::cos(mapFLowToInt[N] * Phi);
                    Q_y_C += pT * std::sin(mapFLowToInt[N] * Phi);
                }
            }
        }
        params_vec.Q_x[N]      = (Q_x);
        params_vec.Q_y[N]      = (Q_y);
        params_vec.Psi_1[N]    = (Psi(Q_y, Q_x) / mapFLowToInt[N]);
        params_vec.Q_x_A[N]    = (Q_x_A);
        params_vec.Q_y_A[N]    = (Q_y_A);
        params_vec.Q_x_B[N]    = (Q_x_B);
        params_vec.Q_y_B[N]    = (Q_y_B);
        params_vec.Q_x_C[N]    = (Q_x_C);
        params_vec.Q_y_C[N]    = (Q_y_C);
        const double Phi_1     = Psi(Q_y_A, Q_x_A) / mapFLowToInt[N];
        const double Phi_CD    = Psi(Q_y_B, Q_x_B) / mapFLowToInt[N];
        const double Phi_TPC   = Psi(Q_y_C, Q_x_C) / mapFLowToInt[N];
        params_vec.cos_N_AB[N] = (std::cos(mapFLowToInt[N] * (Phi_1 - Phi_CD)));
        params_vec.cos_N_AC[N] = (std::cos(mapFLowToInt[N] * (Phi_1 - Phi_TPC)));
        params_vec.cos_N_BC[N] = (std::cos(mapFLowToInt[N] * (Phi_CD - Phi_TPC)));
    }

    return params_vec;
}

auto Coal::getResolutionMap_vec(const EventsMap                &allEvents,
                                const std::vector<std::string> &N_list) -> Coal::ResParamsMap_vec
{
    ResParamsMap_vec ResParamsMap_vec;

    std::map<std::string, double> cos_N_AB;
    std::map<std::string, double> cos_N_AC;
    std::map<std::string, double> cos_N_BC;

    for (const auto &[eventId, OneEvent]: allEvents) {
        ResParamsMap_vec.ResMap[eventId]        = getResolution_vec(OneEvent, N_list);
        ResParamsMap_vec.eventPlaneMap[eventId] = ResParamsMap_vec.ResMap[eventId].Psi_1;
        for (const auto &N: N_list) {
            cos_N_AB[N] += ResParamsMap_vec.ResMap[eventId].cos_N_AB[N];
            cos_N_AC[N] += ResParamsMap_vec.ResMap[eventId].cos_N_AC[N];
            cos_N_BC[N] += ResParamsMap_vec.ResMap[eventId].cos_N_BC[N];
        }
    }
    for (const auto &N: N_list) {
        ResParamsMap_vec.Ref_vn[N] = 1.0;
        if (cos_N_AB[N] * cos_N_AC[N] / cos_N_BC[N] > 0) {
            ResParamsMap_vec.Ref_vn[N] =
                    std::sqrt(cos_N_AB[N] * cos_N_AC[N] / cos_N_BC[N] / allEvents.size());
        }
        std::cout << "Ref_v" << N << ":" << ResParamsMap_vec.Ref_vn[N] << std::endl;
    }

    return ResParamsMap_vec;
}

auto Coal::getResolutionMap(const EventsMap &allEvents, const int N) -> Coal::ResParamsMap
{
    ResParamsMap ResParamsMap;
    double       cos_N_AB     = 0.0;
    double       cos_N_AC     = 0.0;
    double       cos_N_BC     = 0.0;
    const int    total_events = static_cast<int>(allEvents.size());
    for (const auto &[eventId, OneEvent]: allEvents) {
        ResParamsMap.ResMap[eventId] = getResolution(OneEvent, N);
        cos_N_AB += ResParamsMap.ResMap[eventId].cos_N_AB;
        cos_N_AC += ResParamsMap.ResMap[eventId].cos_N_AC;
        cos_N_BC += ResParamsMap.ResMap[eventId].cos_N_BC;
        ResParamsMap.eventPlaneMap[eventId] = ResParamsMap.ResMap[eventId].Psi_1;
    }
    ResParamsMap.Ref_v1 = 1.0;
    if (cos_N_AB * cos_N_AC / cos_N_BC > 0) {
        ResParamsMap.Ref_v2 = std::sqrt(cos_N_AB * cos_N_AC / cos_N_BC / total_events);
    }
    else {
        ResParamsMap.Ref_v2 = 1.0;
    }
    std::cout << "toal_events:" << total_events << ", "
              << "cos_N_AB:" << cos_N_AB << ", "
              << "cos_N_AC:" << cos_N_AC << ", "
              << "cos_N_BC:" << cos_N_BC << "\n";
    std::cout << "Ref_v1:" << ResParamsMap.Ref_v1 << ", "
              << "Ref_v2:" << ResParamsMap.Ref_v2 << std::endl;
    return ResParamsMap;
}

void Coal::rotateEventPlane(ParticleTypeMap &OneEvent, const double angle)
{

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
void Coal::roateAllEventPlane(EventsMap &allEvents, const EventPlaneMap &eventPlaneMap)
{
    for (auto &[eventId, OneEvent]: allEvents) {
        rotateEventPlane(OneEvent, -eventPlaneMap.at(eventId));
    }
}

void Coal::reateAllEventPlaneFlow(EventsMap &allEvents, const EventPlaneMapArray &eventPlaneMap)
{
    for (auto &[eventId, OneEvent]: allEvents) {
        rotateEventPlane(OneEvent, -eventPlaneMap.at(eventId).at("v2"));
    }
}

void Coal::outputAllEventPlaneFlow(const EventsMap                &allEvents,
                                   const EventPlaneMapArray       &eventPlaneMap,
                                   const std::vector<std::string> &FlowName,
                                   const std::string              &eventPlaneName)
{
    std::ofstream output(eventPlaneName);
    for (const auto &eventId: allEvents | std::views::keys) {
        output << eventId << " ";
        for (const auto &name: FlowName) {
            output << eventPlaneMap.at(eventId).at(name) << " ";
        }
        output << "\n";
    }
    output.close();
}

void Coal::formulateFlow(flowResult &result, const double rapidity, const double pT,
                         const double Psi, const double Phi, const RapidityArray &rapidityArray,
                         const std::pair<double, int> &ptBins, const double weight,
                         const std::pair<double, double> &ptRange, const int total_events)
{
    for (const auto &rap: rapidityArray) {
        if (rap.first < rapidity && rapidity <= rap.second) {
            if (const int index = static_cast<int>(pT / ptBins.first);
                index < ptBins.second && index >= 0) {
                result.ptArray[rap][index] += weight / total_events;
                result.v1Array[rap][index] += std::cos(Phi - Psi);
                result.v1Count[rap][index]++;
                result.v2Array[rap][index] += std::cos(2 * (Phi - Psi));
                result.v2Count[rap][index]++;
            }
            if (pT >= ptRange.first && pT <= ptRange.second) {
                result.v1overpt[rap] += std::cos(Phi - Psi);
                result.v1overptCount[rap]++;
                result.v2overpt[rap] += std::cos(2 * (Phi - Psi));
                result.v2overptCount[rap]++;
            }
            result.yeildArray[rap] += weight / total_events;
        }
    }
}

void Coal::formulateFlow_vec(flowResult_vec &result, const double rapidity, const double pT,
                             const std::map<std::string, double> &Psi,
                             const std::vector<std::string> &FlowName, const double Phi,
                             const RapidityArray          &rapidityArray,
                             const std::pair<double, int> &ptBins, const double weight,
                             const std::pair<double, double> &ptRange, const int total_events)
{

    for (const auto &rap: rapidityArray) {
        if (rap.first < rapidity && rapidity <= rap.second) {
            if (const int index = static_cast<int>(pT / ptBins.first);
                index < ptBins.second && index >= 0) {
                result.ptArray[rap][index] += weight / total_events;
                for (const auto &name: FlowName) {
                    result.flowArray[name][rap][index] +=
                            std::cos(mapFLowToInt[name] * (Phi - Psi.at(name)));
                    result.flowCount[name][rap][index]++;
                }
            }
            if (pT >= ptRange.first && pT <= ptRange.second) {
                for (const auto &name: FlowName) {
                    result.flowOverPt[name][rap] +=
                            std::cos(mapFLowToInt[name] * (Phi - Psi.at(name)));
                    result.flowOverPtCount[name][rap]++;
                }
            }
            result.yeildArray[rap] += weight / total_events;
        }
    }
}

auto Coal::calculateFlow(const EventsMap &allEvents, const int pdg,
                         const RapidityArray &rapidityArray, const std::pair<double, int> &ptBins,
                         const ResParamsMap &resolution) -> Coal::flowResult
{

    flowResult          result(rapidityArray, ptBins.second);
    const int           total_events = static_cast<int>(allEvents.size());
    constexpr std::pair ptRange      = {0.4, 2.0};

    result.ptRange = ptRange;
    result.ptBins  = ptBins;

    for (const auto &[eventId, OneEvent]: allEvents) {
        for (const auto &[particleType, particles]: OneEvent) {
            if (particleType == pdg) {
                for (const auto &particle: particles) {
                    const double rapidity = particle.getPseudoRapidity();
                    const double pT       = particle.getPT();
                    const double Phi      = Psi(particle.py, particle.px);
                    const double Psi      = resolution.eventPlaneMap.at(eventId);

                    formulateFlow(result, rapidity, pT, Psi, Phi, rapidityArray, ptBins, 1.0,
                                  ptRange, total_events);
                }
            }
        }
    }
    return result;
}

auto Coal::calculateFlow_vec(const EventsMap &allEvents, int pdg,
                             const RapidityArray            &rapidityArray,
                             const std::pair<double, int>   &ptBins,
                             const std::vector<std::string> &FlowName,
                             const ResParamsMap_vec         &resolution) -> Coal::flowResult_vec
{
    // std::vector<int> FlowName = {1, 2, 3, 4};

    flowResult_vec result(rapidityArray, ptBins.second, FlowName);

    const int total_events = static_cast<int>(allEvents.size());

    constexpr std::pair ptRange = {0.4, 2.0};


    result.ptRange = ptRange;

    result.ptBins = ptBins;


    for (const auto &[eventId, OneEvent]: allEvents) {

        for (const auto &[particleType, particles]: OneEvent) {

            if (particleType == pdg) {

                for (const auto &particle: particles) {


                    const double rapidity = particle.getPseudoRapidity();

                    const double pT = particle.getPT();


                    const double Phi = Psi(particle.py, particle.px);


                    auto Psi = resolution.eventPlaneMap.at(eventId);


                    formulateFlow_vec(result, rapidity, pT, Psi, FlowName, Phi, rapidityArray,
                                      ptBins,

                                      1.0, ptRange, total_events);
                }
            }
        }
    }

    return result;
}


auto Coal::calculateFlowMatrix(const Eigen::MatrixXd &targetParticles,
                               const RapidityArray &rapidityArray, const ClusterParams &params,
                               const ResParamsMap &resolution) -> Coal::flowResult
{
    const int                 ptBins      = params.ptBins.second;
    const int                 totalEvents = static_cast<int>(params.eventFactor);
    std::pair<double, double> ptRange;
    if (params.NBody == 2) {
        ptRange = {0.8, 2.0};
    }
    else if (params.NBody == 3) {
        ptRange = {1.2, 3.0};
    }
    else if (params.NBody == 4) {
        ptRange = {1.6, 4.0};
    }
    else if (params.NBody == 8) {
        ptRange = {3.2, 8.0};
    }
    else if (params.NBody == 12) {
        ptRange = {4.8, 12.0};
    }
    flowResult result(rapidityArray, ptBins);
    result.ptRange = ptRange;
    result.ptBins  = params.ptBins;
    double psi;

    int particleCount = 0;
    int currentIndex  = 0;

    for (int i = 0; i < targetParticles.rows(); ++i) {
        const double rapidity = getMatrixPseudoRapidity(targetParticles, i);
        const double pT       = getMatrixPt(targetParticles, i);
        const double phi      = Psi(targetParticles(i, 2), targetParticles(i, 1));
        if (params.mixEvents == 1) {
            if (particleCount >= resolution.selectEventID.at(currentIndex).second) {
                currentIndex++;
                particleCount = 0;
            }
            int eventId = resolution.selectEventID.at(currentIndex).first;
            psi         = resolution.eventPlaneMap.at(eventId);
            particleCount++;
        }
        else {
            psi = 0.0;
        }
        formulateFlow(result, rapidity, pT, psi, phi, rapidityArray, params.ptBins,
                      targetParticles(i, 10), ptRange, totalEvents);
    }
    return result;
}

auto Coal::calculateFlowMatrix_vec(const Eigen::MatrixXd &targetParticles,
                                   const RapidityArray &rapidityArray, const ClusterParams &params,
                                   const ResParamsMap_vec         &resolution,
                                   const std::vector<std::string> &FlowName) -> Coal::flowResult_vec
{
    const int ptBins = params.ptBins.second;

    const int totalEvents = static_cast<int>(params.eventFactor);

    std::pair<double, double> ptRange;

    if (params.NBody == 2) {

        ptRange = {0.8, 2.0};
    }
    else if (params.NBody == 3) {

        ptRange = {1.2, 3.0};
    }
    else if (params.NBody == 4) {

        ptRange = {1.6, 4.0};
    }
    else if (params.NBody == 8) {

        ptRange = {3.2, 8.0};
    }
    else if (params.NBody == 12) {

        ptRange = {4.8, 12.0};
    }
    // std::vector<int> FlowName = {1, 2, 3, 4};

    flowResult_vec result(rapidityArray, ptBins, FlowName);
    result.ptRange = ptRange;

    result.ptBins = params.ptBins;

    std::map<std::string, double> psi;
    int                           particleCount = 0;

    int currentIndex = 0;
    for (int i = 0; i < targetParticles.rows(); ++i) {

        const double rapidity = getMatrixPseudoRapidity(targetParticles, i);

        const double pT = getMatrixPt(targetParticles, i);

        const double phi = Psi(targetParticles(i, 2), targetParticles(i, 1));

        if (params.mixEvents == 1) {

            if (particleCount >= resolution.selectEventID.at(currentIndex).second) {

                currentIndex++;

                particleCount = 0;
            }

            int eventId = resolution.selectEventID.at(currentIndex).first;

            psi = resolution.eventPlaneMap.at(eventId);

            particleCount++;
        }
        else {
            for (const auto &name: FlowName) {
                psi[name] = 0.0;
            }
        }

        formulateFlow_vec(result, rapidity, pT, psi, FlowName, phi, rapidityArray, params.ptBins,

                          targetParticles(i, 10), ptRange, totalEvents);
    }

    return result;
}


auto Coal::getMatrixPt(const Eigen::MatrixXd &targetParticles, const int i) -> double
{
    return std::sqrt(targetParticles(i, 1) * targetParticles(i, 1) +
                     targetParticles(i, 2) * targetParticles(i, 2));
}

auto Coal::getMatrixRapidity(const Eigen::MatrixXd &targetParticles, const int i) -> double
{
    const double rapidity = 0.5 * log((targetParticles(i, 4) + targetParticles(i, 3)) /
                                      (targetParticles(i, 4) - targetParticles(i, 3)));
    return rapidity;
}
auto Coal::getMatrixPseudoRapidity(const Eigen::MatrixXd &targetParticles, const int i) -> double
{
    const double p_total = std::sqrt(targetParticles(i, 1) * targetParticles(i, 1) +
                                     targetParticles(i, 2) * targetParticles(i, 2) +
                                     targetParticles(i, 3) * targetParticles(i, 3));

    const double pseudoRapidity =
            0.5 * log((p_total + targetParticles(i, 3)) / (p_total - targetParticles(i, 3)));

    return pseudoRapidity;
}
