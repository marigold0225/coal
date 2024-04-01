//
// Created by mafu on 24-3-26.
//
#include "../headers/flow.h"
#include "gsl/gsl_sf_bessel.h"
#include <ranges>
#include <fstream>
#include <iostream>

int Coal::eventPlaneEquation(const gsl_vector *x, void *p, gsl_vector *f) {
    const params *params = static_cast<struct params *>(p);
    const double a     = params->a;
    const int k        = params->k;
    const double x_val = gsl_vector_get(x, 0);
    const double arg   = x_val * x_val / 4;

    const double I_k_minus_1_over_2 = gsl_sf_bessel_Inu((k - 1) / 2.0, arg);
    const double I_k_plus_1_over_2  = gsl_sf_bessel_Inu((k + 1) / 2.0, arg);

    const double y = M_PI / 2 / sqrt(2) * x_val * exp(-arg) * (I_k_minus_1_over_2 + I_k_plus_1_over_2) - a;

    gsl_vector_set(f, 0, y);
    return GSL_SUCCESS;
}
double Coal::solveEquation(const double a,const int k) {
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *s            = gsl_multiroot_fsolver_alloc(T, 1);

    params par                  = {a, k};
    gsl_multiroot_function func = {&eventPlaneEquation, 1, &par};

     constexpr double x_init[1] = {0.5};
    gsl_vector *x    = gsl_vector_alloc(1);
    gsl_vector_set(x, 0, x_init[0]);

    gsl_multiroot_fsolver_set(s, &func, x);

    int status, iter = 0;
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        if (status)
            break;

        status = gsl_multiroot_test_residual(s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 1000);

    const double root = gsl_vector_get(s->x, 0);
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return root;
}
double Coal::getRef(const double x, const int k) {
    const double arg                = x * x / 4;
    const double I_k_minus_1_over_2 = gsl_sf_bessel_Inu((k - 1) / 2.0, arg);
    const double I_k_plus_1_over_2  = gsl_sf_bessel_Inu((k + 1) / 2.0, arg);
    return M_PI / 2 / sqrt(2) * x * std::exp(-arg) *
           (I_k_minus_1_over_2 + I_k_plus_1_over_2);
}
double Coal::Psi(double y, double x) {
    double angle = atan2(y, x);
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    return angle;
}
Coal::ResParams Coal::getResolution(const ParticleEventMap &OneEvent, int N) {

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
void Coal::rotateEventPlane(ParticleTypeMap &OneEvent, double angle) {
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
void Coal::getFlow(const EventsMap &allEvents, int pdg,
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
