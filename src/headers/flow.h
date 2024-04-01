//
// Created by mafu on 24-3-26.
//
#pragma once
#include "gsl/gsl_multiroots.h"
#include "fileloader.h"

namespace Coal {

    struct params {
        double a;
        int k;
    };

    struct ResParams {
        double Q_x;
        double Q_y;
        double Psi_1;
        double Q_x_A;
        double Q_y_A;
        double Q_x_B;
        double Q_y_B;
        double Q_x_C;
        double Q_y_C;
        double cos_N_AB;
        double cos_N_AC;
        double cos_N_BC;
    };

    int eventPlaneEquation(const gsl_vector *x, void *p, gsl_vector *f);

    double solveEquation(double a, int k);

    double getRef(double x, int k);

    double Psi(double y, double x);

    ResParams getResolution(const ParticleEventMap &OneEvent, int N);

    void rotateEventPlane(ParticleTypeMap &OneEvent, double angle);

    void getFlow(const EventsMap &allEvents, int pdg,
                          const std::string &outputFilename,
                          const RapidityArray &rapidityArray,
                          const std::pair<double, int> &ptBins);

}// namespace Coal
