//
// Created by mafu on 24-3-26.
//
#pragma once
#include "config.h"
#include "shortcut.h"
#include <vector>
// #include "gsl/gsl_multiroots.h"
// #include "gsl/gsl_sf_bessel.h"
namespace Coal {

    struct params
    {
        double a;
        int    k;
    };

    struct ResParams
    {
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
    inline std::map<std::string, int> mapFLowToInt = {
            {"v1", 1},
            {"v2", 2},
            {"v3", 3},
            {"v4", 4},
    };
    struct ResParams_vec
    {

        FlowType Q_x{};

        FlowType Q_y{};

        FlowType Psi_1{};

        FlowType Q_x_A{};

        FlowType Q_y_A{};

        FlowType Q_x_B{};

        FlowType Q_y_B{};

        FlowType Q_x_C{};

        FlowType Q_y_C{};

        FlowType cos_N_AB{};

        FlowType cos_N_AC{};

        FlowType cos_N_BC{};
    };
    struct ResParamsMap_vec
    {

        std::map<int, ResParams_vec> ResMap;

        EventPlaneMapArray eventPlaneMap;

        FlowType Ref_vn{};

        std::map<int, Pair> selectEventID;
    };
    struct ResParamsMap
    {
        std::map<int, ResParams> ResMap;
        EventPlaneMap            eventPlaneMap;
        double                   Ref_v1;
        double                   Ref_v2;
        std::map<int, Pair>      selectEventID;
    };

    struct flowResult
    {
        std::map<RapidityRange, std::vector<double>> ptArray;
        std::map<RapidityRange, double>              yeildArray;
        std::map<RapidityRange, std::vector<double>> v1Array;
        std::map<RapidityRange, std::vector<double>> v2Array;
        std::map<RapidityRange, std::vector<int>>    v1Count;
        std::map<RapidityRange, std::vector<int>>    v2Count;
        std::map<RapidityRange, double>              v1overpt;
        std::map<RapidityRange, int>                 v1overptCount;
        std::map<RapidityRange, double>              v2overpt;
        std::map<RapidityRange, int>                 v2overptCount;
        std::pair<double, double>                    ptRange;
        std::pair<double, int>                       ptBins;

        explicit flowResult(const RapidityArray &rapidityArray, int bins);
    };

    template<typename T> using FlowMap = std::map<RapidityRange, std::vector<T>>;

    struct flowResult_vec
    {

        std::map<RapidityRange, std::vector<double>> ptArray;

        std::map<RapidityRange, double> yeildArray;

        std::map<std::string, FlowMap<double>> flowArray;

        std::map<std::string, FlowMap<int>> flowCount;

        std::map<std::string, std::map<RapidityRange, double>> flowOverPt;

        std::map<std::string, std::map<RapidityRange, int>> flowOverPtCount;

        std::pair<double, double> ptRange;

        std::pair<double, int> ptBins;


        explicit flowResult_vec(const RapidityArray &rapidityArray, int bins,
                                const std::vector<std::string> &FlowName);
    };

    int eventPlaneEquation_std(double x_val, const params *p, double &y);

    double solveEquation_std(double a, int k);

    double getRef_std(double x, int k);
    double Psi(double y, double x);

    ResParams     getResolution(const ParticleEventMap &OneEvent, int N);
    ResParams_vec getResolution_vec(const ParticleEventMap         &OneEvent,
                                    const std::vector<std::string> &N_list);

    ResParamsMap     getResolutionMap(const EventsMap &allEvents, int N);
    ResParamsMap_vec getResolutionMap_vec(const EventsMap                &allEvents,
                                          const std::vector<std::string> &N_list);

    void rotateEventPlane(ParticleTypeMap &OneEvent, double angle);

    void roateAllEventPlane(EventsMap &allEvents, const EventPlaneMap &eventPlaneMap);

    void reateAllEventPlaneFlow(EventsMap &allEvents, const EventPlaneMapArray &eventPlaneMap);

    void outputAllEventPlaneFlow(const EventsMap          &allEvents,
                                 const EventPlaneMapArray &eventPlaneMap,
                                 const std::vector<std::string> &FlowName,
                                 const std::string & eventPlaneName);

    void formulateFlow(flowResult &result, double rapidity, double pT, double Psi, double Phi,
                       const RapidityArray &rapidityArray, const std::pair<double, int> &ptBins,
                       double weight, const std::pair<double, double> &ptRange, int total_events);


    void formulateFlow_vec(flowResult_vec &result, double rapidity, double pT,
                           const std::map<std::string, double> &Psi,
                           const std::vector<std::string>      &FlowName,

                           double Phi, const RapidityArray &rapidityArray,

                           const std::pair<double, int> &ptBins, double weight,

                           const std::pair<double, double> &ptRange, int total_events);

    flowResult calculateFlow(const EventsMap &allEvents, int pdg,
                             const RapidityArray          &rapidityArray,
                             const std::pair<double, int> &ptBins, const ResParamsMap &resolution);


    flowResult_vec calculateFlow_vec(const EventsMap &allEvents, int pdg,

                                     const RapidityArray &rapidityArray,

                                     const std::pair<double, int>   &ptBins,
                                     const std::vector<std::string> &FlowName,
                                     const ResParamsMap_vec         &resolution);


    flowResult calculateFlowMatrix(const Eigen::MatrixXd &targetParticles,
                                   const RapidityArray &rapidityArray, const ClusterParams &params,
                                   const ResParamsMap &resolution);

    flowResult_vec calculateFlowMatrix_vec(const Eigen::MatrixXd &targetParticles,

                                           const RapidityArray &rapidityArray,

                                           const ClusterParams &params,

                                           const ResParamsMap_vec         &resolution,
                                           const std::vector<std::string> &FlowName);

    static double getMatrixPt(const Eigen::MatrixXd &targetParticles, int i);

    static double getMatrixRapidity(const Eigen::MatrixXd &targetParticles, int i);

    static double getMatrixPseudoRapidity(const Eigen::MatrixXd &targetParticles, int i);


}// namespace Coal
