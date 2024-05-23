//
// Created by mafu on 4/2/2024.
//
#include "../headers/config.h"

#include <iostream>
#include <numeric>

Coal::ConfigParser::ConfigParser(const std::string &filename) {
    YAML::Node config = YAML::LoadFile(filename);
    total             = config;
    mode              = config["General"]["Mode"].as<std::string>();
    general           = config["General"];
    inputFileName     = config["General"]["Path"].as<std::string>() + "/" +
                    config["General"]["Filename"].as<std::string>();
    rapidityArray = config["Output"]["RapidityRange"].as<RapidityArray>();
    output        = config["Output"];
    seed          = config["General"]["Seed"].as<int>();
    outputPath    = config["Output"]["Path"].as<std::string>();
    clusterParams = config["ClusterParams"];

    for (const auto &reaction: config["Reactions"]) {
        reactionSwitch[reaction.first.as<std::string>()] = reaction.second.as<bool>();
    }
}

Coal::ClusterParams::ClusterParams(const YAML::Node &node) {
    NBody       = node["NBody"].as<int>();
    Loop        = node["Loop"].as<int>();
    mixEvents   = node["MixEvents"].as<int>();
    eventFactor = Loop * pow(mixEvents, NBody);
    probFactor  = pow(8, NBody - 1);
    precision   = node["Precision"].as<int>();
    if (const auto gc_value = node["gc"].as<std::string>();
        gc_value.find('/') != std::string::npos) {
        gc = parseFraction(gc_value);
    } else {
        gc = node["gc"].as<double>();
    }
    pdg            = node["PDG"].as<int>();
    probabilityCut = node["ProbCut"].as<double>();
    Fixed          = node["MassArray"]["Fixed"].as<bool>();
    if (Fixed) {
        MassArray = node["MassArray"]["Array"].as<std::vector<double>>();
        setMassArray(MassArray);
    } else {
        MassArray = node["MassArray"]["Array"].as<std::vector<double>>();
    }
    // MassArray      = node["MassArray"].as<std::vector<double>>();
    SigArray       = node["Sig"].as<std::vector<double>>();
    PDGArray       = node["From"].as<std::vector<int>>();
    ptBins         = node["Pt"].as<std::pair<double, int>>();
    originRapidity = node["RapidityCut"].as<RapidityArray>();
    targetRapidity = node["TargetRapidity"].as<RapidityRange>();
    setMassArray(MassArray);
}

void Coal::ClusterParams::setMassArray(const std::vector<double> &massArray) {
    M.resize(massArray.size() - 1);
    M_inv_t.resize(massArray.size() - 1);
    for (int i = 1; i < massArray.size(); ++i) {
        std::vector temMassArray(massArray.begin(), massArray.begin() + i + 1);
        M[i - 1]       = MassMatrix(temMassArray);
        M_inv_t[i - 1] = M[i - 1].inverse().transpose();
    }
}

double Coal::parseFraction(const std::string &fraction) {
    std::istringstream iss(fraction);
    double numerator, denominator;
    if (char slash; !(iss >> numerator >> slash >> denominator)) {
        std::cerr << "Error: cannot parse fraction " << fraction << std::endl;
        exit(1);
    }
    return numerator / denominator;
}

Eigen::MatrixXd Coal::MassMatrix(const std::vector<double> &massArray) {
    const int n       = static_cast<int>(massArray.size());
    // Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);

    const double totalMass = Eigen::VectorXd::Map(massArray.data(), n).sum();

    M.row(0) = Eigen::VectorXd::Map(massArray.data(), n) / totalMass;

    for (int k = 1; k < n; ++k) {
        if (k == 1) {
            M(k, 0) = 1.0 / sqrt(2);
            M(k, 1) = -1.0 / sqrt(2);
        } else {
            double sum = 0;
            for (int i = 0; i < k; ++i) {
                sum += massArray[i];
            }
            for (int i = 0; i < k; ++i) {
                M(k, i) = sqrt(static_cast<double>(k) / (k + 1)) * massArray[i] / sum;
            }
            M(k, k) = -sqrt(static_cast<double>(k) / (k + 1));
        }
    }

    return M;
}