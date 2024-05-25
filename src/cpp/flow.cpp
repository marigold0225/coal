//
// Created by mafu on 24-3-26.
//
#include "../headers/flow.h"
#include <fstream>
#include <iostream>
#include <ranges>

Coal::flowResult::flowResult(const RapidityArray &rapidityArray, const int bins) {
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

// int Coal::eventPlaneEquation(const gsl_vector *x, void *p, gsl_vector *f) {
//     const params *params = static_cast<struct params *>(p);
//     const double a       = params->a;
//     const int k          = params->k;
//     const double x_val   = gsl_vector_get(x, 0);
//     const double arg     = x_val * x_val / 4;
//
//     const double I_k_minus_1_over_2 = gsl_sf_bessel_Inu((k - 1) / 2.0, arg);
//     // const double I = std::cyl_bessel_i()
//     const double I_k_plus_1_over_2 = gsl_sf_bessel_Inu((k + 1) / 2.0, arg);
//
//     const double y = M_PI / 2 / sqrt(2) * x_val * exp(-arg) *
//                              (I_k_minus_1_over_2 + I_k_plus_1_over_2) -
//                      a;
//
//     gsl_vector_set(f, 0, y);
//     return GSL_SUCCESS;
// }
int Coal::eventPlaneEquation_std(const double x_val, const params *p, double &y) {
    const double a   = p->a;
    const int k      = p->k;
    const double arg = x_val * x_val / 4;

    const double I_k_minus_1_over_2 = std::cyl_bessel_i((k - 1) / 2.0, arg);
    const double I_k_plus_1_over_2  = std::cyl_bessel_i((k + 1) / 2.0, arg);

    y = M_PI / 2 / sqrt(2) * x_val * exp(-arg) *
                (I_k_minus_1_over_2 + I_k_plus_1_over_2) -
        a;
    return 0;
}
// double Coal::solveEquation(const double a, const int k) {
//     const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
//     gsl_multiroot_fsolver *s            = gsl_multiroot_fsolver_alloc(T, 1);
//
//     params par                  = {a, k};
//     gsl_multiroot_function func = {&eventPlaneEquation, 1, &par};
//
//     constexpr double x_init[1] = {0.5};
//     gsl_vector *x              = gsl_vector_alloc(1);
//     gsl_vector_set(x, 0, x_init[0]);
//
//     gsl_multiroot_fsolver_set(s, &func, x);
//
//     int status, iter = 0;
//     do {
//         iter++;
//         status = gsl_multiroot_fsolver_iterate(s);
//         if (status)
//             break;
//
//         status = gsl_multiroot_test_residual(s->f, 1e-7);
//     } while (status == GSL_CONTINUE && iter < 1000);
//
//     const double root = gsl_vector_get(s->x, 0);
//     gsl_multiroot_fsolver_free(s);
//     gsl_vector_free(x);
//
//     return root;
// }
double Coal::solveEquation_std(const double a, const int k) {
    const params par       = {a, k};
    double x               = 0.5;
    constexpr int max_iter = 1000;
    for (int iter = 0; iter < max_iter; ++iter) {
        double y;
        eventPlaneEquation_std(x, &par, y);
        if (constexpr double tol = 1e-7; std::abs(y) < tol)
            break;
        constexpr double h = 1e-5;
        double y1, y2;
        eventPlaneEquation_std(x + h, &par, y1);
        eventPlaneEquation_std(x - h, &par, y2);
        const double dy = (y1 - y2) / (2 * h);
        x -= y / dy;
    }
    return x;
}
// double Coal::getRef(const double x, const int k) {
//     const double arg                = x * x / 4;
//     const double I_k_minus_1_over_2 = gsl_sf_bessel_Inu((k - 1) / 2.0, arg);
//     const double I_k_plus_1_over_2  = gsl_sf_bessel_Inu((k + 1) / 2.0, arg);
//     return M_PI / 2 / sqrt(2) * x * std::exp(-arg) *
//            (I_k_minus_1_over_2 + I_k_plus_1_over_2);
// }
double Coal::getRef_std(const double x, const int k) {
    const double arg                = x * x / 4;
    const double I_k_minus_1_over_2 = std::cyl_bessel_i((k - 1) / 2.0, arg);
    const double I_k_plus_1_over_2  = std::cyl_bessel_i((k + 1) / 2.0, arg);
    return M_PI / 2 / sqrt(2) * x * std::exp(-arg) *
           (I_k_minus_1_over_2 + I_k_plus_1_over_2);
}

double Coal::Psi(const double y, const double x) {
    double angle = atan2(y, x);
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    return angle;
}

Coal::ResParams Coal::getResolution(const ParticleEventMap &OneEvent, const int N) {

    ResParams params{};
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
            // const double rapidity       = particle.getRapidity();
            const double pT  = particle.getPT();
            // const double Phi = particle.getPhi();
            const double Phi = Psi(particle.py, particle.px);
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
            if (pseudoRapidity > -0.8 && pseudoRapidity < 0.8 && pT > 0.2 && pT < 20) {
                Q_x += pT * std::cos(2 * Phi);
                Q_y += pT * std::sin(2 * Phi);
            }
            if (pseudoRapidity > -3.7 && pseudoRapidity < -1.7 ||
                pseudoRapidity > 2.8 && pseudoRapidity < 5.1) {
                Q_x_A += pT * std::cos(2 * Phi);
                Q_y_A += pT * std::sin(2 * Phi);
            }
            if (pseudoRapidity < 0.8 && pseudoRapidity > 0.0 && pT > 0.2 && pT < 20) {
                Q_x_B += pT * std::cos(2 * Phi);
                Q_y_B += pT * std::sin(2 * Phi);
            }
            if (pseudoRapidity > -0.8 && pseudoRapidity < 0.0 && pT > 0.2 && pT < 20) {
                Q_x_C += pT * std::cos(2 * Phi);
                Q_y_C += pT * std::sin(2 * Phi);
            }
        }
    }
    params.Q_x           = Q_x;
    params.Q_y           = Q_y;
    params.Psi_1         = Psi(Q_y, Q_x);
    params.Q_x_A         = Q_x_A;
    params.Q_y_A         = Q_y_A;
    params.Q_x_B         = Q_x_B;
    params.Q_y_B         = Q_y_B;
    params.Q_x_C         = Q_x_C;
    params.Q_y_C         = Q_y_C;
    const double Phi_1   = Psi(Q_y_A, Q_x_A);
    const double Phi_CD  = Psi(Q_y_B, Q_x_B);
    const double Phi_TPC = Psi(Q_y_C, Q_x_C);

    params.cos_N_AB = std::cos(N * (Phi_1 - Phi_CD));
    params.cos_N_AC = std::cos(N * (Phi_1 - Phi_TPC));
    params.cos_N_BC = std::cos(N * (Phi_CD - Phi_TPC));
    return params;
}
Coal::ResParamsMap Coal::getResolutionMap(const EventsMap &allEvents, const int N) {
    ResParamsMap ResParamsMap;
    double cos_N_AB        = 0.0;
    double cos_N_AC        = 0.0;
    double cos_N_BC        = 0.0;
    const int total_events = static_cast<int>(allEvents.size());
    for (const auto &[eventId, OneEvent]: allEvents) {
        ResParamsMap.ResMap[eventId] = getResolution(OneEvent, N);
        cos_N_AB += ResParamsMap.ResMap[eventId].cos_N_AB;
        cos_N_AC += ResParamsMap.ResMap[eventId].cos_N_AC;
        cos_N_BC += ResParamsMap.ResMap[eventId].cos_N_BC;
        ResParamsMap.eventPlaneMap[eventId] = ResParamsMap.ResMap[eventId].Psi_1;
    }
    // ResParamsMap.Ref_v1  = std::sqrt(cos_N_AB * cos_N_AC / cos_N_BC / total_events);
    // ResParamsMap.Ref_v2 = std::sqrt(cos_N_AB * cos_N_AC / cos_N_BC / total_events);
    ResParamsMap.Ref_v1 = 0.0;
    if (cos_N_AB * cos_N_AC / cos_N_BC > 0) {
        ResParamsMap.Ref_v2 = std::sqrt(cos_N_AB * cos_N_AC / cos_N_BC / total_events);
    } else {
        ResParamsMap.Ref_v2 = 0.0;
    }
    std::cout << "toal_events:" << total_events << ", "
              << "cos_N_AB:" << cos_N_AB << ", "
              << "cos_N_AC:" << cos_N_AC << ", "
              << "cos_N_BC:" << cos_N_BC << "\n";
    std::cout << "Ref_v1:" << ResParamsMap.Ref_v1 << ", "
              << "Ref_v2:" << ResParamsMap.Ref_v2 << std::endl;
    // ResParamsMap.Ref_v2  = getRef(solveEquation(ResParamsMap.Ref_v1, 1), 2);
    return ResParamsMap;
}

void Coal::rotateEventPlane(ParticleTypeMap &OneEvent, const double angle) {

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
void Coal::roateAllEventPlane(EventsMap &allEvents, const EventPlaneMap &eventPlaneMap) {
    for (auto &[eventId, OneEvent]: allEvents) {
        rotateEventPlane(OneEvent, -eventPlaneMap.at(eventId));
    }
}
void Coal::formulateFlow(flowResult &result, const double rapidity, const double pT,
                         const double Psi, const double Phi,
                         const RapidityArray &rapidityArray,
                         const std::pair<double, int> &ptBins, const double weight,
                         const std::pair<double, double> &ptRange,
                         const int total_events) {
    for (auto &rap: rapidityArray) {
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

Coal::flowResult Coal::calculateFlow(const EventsMap &allEvents, const int pdg,
                                     const RapidityArray &rapidityArray,
                                     const std::pair<double, int> &ptBins,
                                     const ResParamsMap &resolution) {

    flowResult result(rapidityArray, ptBins.second);
    const int total_events      = static_cast<int>(allEvents.size());
    constexpr std::pair ptRange = {0.4, 2.0};

    result.ptRange = ptRange;
    result.ptBins  = ptBins;

    for (const auto &[eventId, OneEvent]: allEvents) {
        for (const auto &[particleType, particles]: OneEvent) {
            if (particleType == pdg) {
                for (const auto &particle: particles) {
                    const double rapidity = particle.getRapidity();
                    // const double rapidity = particle.getPseudoRapidity();
                    const double pT  = particle.getPT();
                    // const double Phi = particle.getPhi();
                    const double Phi = Psi(particle.py,particle.px);
                    // constexpr double Psi  = 0.0;
                    const double Psi = resolution.eventPlaneMap.at(eventId);

                    formulateFlow(result, rapidity, pT, Psi, Phi, rapidityArray, ptBins,
                                  1.0, ptRange, total_events);
                }
            }
        }
    }
    return result;
}
Coal::flowResult Coal::calculateFlowMatrix(const Eigen::MatrixXd &targetParticles,
                                           const RapidityArray &rapidityArray,
                                           const ClusterParams &params,
                                           const ResParamsMap &resolution) {
    const int ptBins      = params.ptBins.second;
    const int totalEvents = static_cast<int>(params.eventFactor);
    std::pair<double, double> ptRange;
    if (params.NBody == 2) {
        ptRange = {0.8, 2.0};
    } else if (params.NBody == 3) {
        ptRange = {1.2, 3.0};
    } else if (params.NBody == 4) {
        ptRange = {1.6, 4.0};
    } else if (params.NBody == 8) {
        ptRange = {3.2, 8.0};
    } else if (params.NBody == 12) {
        ptRange = {4.8, 12.0};
    }
    flowResult result(rapidityArray, ptBins);
    result.ptRange = ptRange;
    result.ptBins  = params.ptBins;
    double psi;

    int particleCount = 0;
    int currentIndex  = 0;

    for (int i = 0; i < targetParticles.rows(); ++i) {
        // const double rapidity = getMatrixRapidity(targetParticles, i);
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
        } else {
            psi = 0.0;
        }
        formulateFlow(result, rapidity, pT, psi, phi, rapidityArray, params.ptBins,
                      targetParticles(i, 10), ptRange, totalEvents);
    }
    return result;
}

double Coal::getMatrixPt(const Eigen::MatrixXd &targetParticles, const int i) {
    return std::sqrt(targetParticles(i, 1) * targetParticles(i, 1) +
                     targetParticles(i, 2) * targetParticles(i, 2));
}

double Coal::getMatrixRapidity(const Eigen::MatrixXd &targetParticles, const int i) {
    const double rapidity = 0.5 * log((targetParticles(i, 4) + targetParticles(i, 3)) /
                                      (targetParticles(i, 4) - targetParticles(i, 3)));
    return rapidity;
}
double Coal::getMatrixPseudoRapidity(const Eigen::MatrixXd &targetParticles,
                                     const int i) {
    const double p_total = std::sqrt(targetParticles(i, 1) * targetParticles(i, 1) +
                                     targetParticles(i, 2) * targetParticles(i, 2) +
                                     targetParticles(i, 3) * targetParticles(i, 3));

    const double pseudoRapidity = 0.5 * log((p_total + targetParticles(i, 3)) /
                                            (p_total - targetParticles(i, 3)));

    return pseudoRapidity;
}
