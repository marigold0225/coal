//
// Created by mafu on 24-3-26.
//
#pragma once
#include "config.h"
// #include "gsl/gsl_multiroots.h"
// #include "gsl/gsl_sf_bessel.h"
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

    struct ResParamsMap {
        std::map<int, ResParams> ResMap;
        EventPlaneMap eventPlaneMap;
        double Ref_v1;
        double Ref_v2;
        std::map<int,Pair> selectEventID;
    };

    struct flowResult {
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
        std::pair<double, double> ptRange;
        std::pair<double, int> ptBins;

        explicit flowResult(const RapidityArray &rapidityArray,int bins);
    };

    // int eventPlaneEquation(const gsl_vector *x, void *p, gsl_vector *f);
    int eventPlaneEquation_std(double x_val, const params *p, double &y);

    // double solveEquation(double a, int k);
    double solveEquation_std(double a, int k);

    // double getRef(double x, int k);
    double getRef_std(double x, int k);
    double Psi(double y, double x);

    ResParams getResolution(const ParticleEventMap &OneEvent, int N);

    ResParamsMap getResolutionMap(const EventsMap &allEvents, int N);

    void rotateEventPlane(ParticleTypeMap &OneEvent, double angle);

    void roateAllEventPlane(EventsMap &allEvents, const EventPlaneMap &eventPlaneMap);

    void formulateFlow(flowResult& result, double rapidity,double pT,double Psi, double Phi,
                        const RapidityArray &rapidityArray, const std::pair<double, int> &ptBins,
                        double weight,const std::pair<double, double> &ptRange,int total_events);

    flowResult calculateFlow(const EventsMap &allEvents, int pdg,
                             const RapidityArray &rapidityArray,
                             const std::pair<double, int> &ptBins,
                             const ResParamsMap &resolution);

    flowResult calculateFlowMatrix(const Eigen::MatrixXd &targetParticles,
                                   const RapidityArray &rapidityArray,
                                   const ClusterParams &params,
                                   const ResParamsMap &resolution);


    static double getMatrixPt(const Eigen::MatrixXd &targetParticles, int i);

    static double getMatrixRapidity(const Eigen::MatrixXd &targetParticles, int i);

    static double getMatrixPseudoRapidity(const Eigen::MatrixXd &targetParticles, int i);


}// namespace Coal
