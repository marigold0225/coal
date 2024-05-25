//
// Created by mafu on 4/2/2024.
//
#pragma once

// #define EIGEN_USE_MKL_ALL
#include "shortcut.h"
#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

namespace Coal {
    struct ConfigParser {
        YAML::Node total;
        std::string mode;
        std::string inputFileName;
        std::string outputPath;
        RapidityArray rapidityArray;
        std::map<std::string, bool> reactionSwitch;
        YAML::Node clusterParams;
        int seed;
        YAML::Node general;
        YAML::Node output;
        explicit ConfigParser(const std::string &filename);
    };

    struct ClusterParams {
        int NBody;
        int Loop;
        int mixEvents;
        double eventFactor;
        double probFactor;
        int precision;
        double gc;
        int pdg;
        bool Fixed;
        std::vector<double> MassArray;
        std::vector<double> SigArray;
        std::vector<double> Sig_inv_sq;
        std::vector<double> Sig_sq_hbar2;
        std::vector<int> PDGArray;
        double probabilityCut;
        RapidityArray originRapidity;
        RapidityRange targetRapidity;
        std::pair<double, int> ptBins;

        std::vector<Eigen::MatrixXd> M;
        std::vector<Eigen::MatrixXd> M_inv_t;
        explicit ClusterParams(const YAML::Node &node);

        void setMassArray(const std::vector<double> &massArray);
    };

    double parseFraction(const std::string &fraction);

    Eigen::MatrixXd MassMatrix(const std::vector<double> &massArray);
}

