//
// Created by mafu on 1/6/2024.
//

#pragma once
#include "smash.h"
#include <Eigen/Dense>
#include <utility>

namespace Coal {

    struct ClusterParams {
        int NBody;
        int Loop;
        int mixEvents;
        int nmix;
        double eventFactor;
        double probFactor;
        int precision;
        double gc;
        int pdg;
        std::vector<double> MassArray;
        std::vector<double> SigArray;
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

    RapidityArray defineRapidityRange(const YAML::Node &node);

    bool checkFileExits(const std::string &filename,
                        const std::vector<std::string> &labels,
                        const std::string &fileType);

    std::string constructFilename(const std::string &filename,
                                  const std::string &fileType,
                                  const std::string &label);
    std::vector<double> linspace(double start, double end, int num);

    EventsMap selectEvents(const EventsMap &eventMap,
                           const ClusterParams &params);

    MultiParticleArray selectParticles(const EventsMap &eventMap,
                                       const ClusterParams &params);
    std::vector<Eigen::MatrixXd>
    selectParticlesMatrix(const EventsMap &eventMap,
                          const ClusterParams &params);

    Eigen::MatrixXd MassMatrix(const std::vector<double> &massArray);

    ParticleArray boostToCOM(const ParticleArray &particles,
                             const Particle &targetParticle);
    //martix boost
    void lorenzBoostMatrix(Eigen::MatrixXd &particles, double beta_x,
                           double beta_y, double beta_z);

    Eigen::MatrixXd boostToComMatrix(const Eigen::MatrixXd &particles);

    std::vector<Eigen::Matrix4d> calculateLorentz(const Eigen::VectorXd &betaX,
                                                  const Eigen::VectorXd &betaY,
                                                  const Eigen::VectorXd &betaZ);

    void applyLorentzBoost(Eigen::Ref<Eigen::MatrixXd> combinedX,
                           Eigen::Ref<Eigen::MatrixXd> combinedY,
                           Eigen::Ref<Eigen::MatrixXd> combinedZ,
                           Eigen::Ref<Eigen::MatrixXd> combinedT,
                           Eigen::Ref<Eigen::MatrixXd> combinedPX,
                           Eigen::Ref<Eigen::MatrixXd> combinedPY,
                           Eigen::Ref<Eigen::MatrixXd> combinedPZ,
                           Eigen::Ref<Eigen::MatrixXd> combinedP0,
                           const std::vector<Eigen::Matrix4d> &lorentz);


    void applyLorentzBoostCols(Eigen::Ref<Eigen::MatrixXd> combinedX,
                               Eigen::Ref<Eigen::MatrixXd> combinedY,
                               Eigen::Ref<Eigen::MatrixXd> combinedZ,
                               Eigen::Ref<Eigen::MatrixXd> combinedT,
                               Eigen::Ref<Eigen::MatrixXd> combinedPX,
                               Eigen::Ref<Eigen::MatrixXd> combinedPY,
                               Eigen::Ref<Eigen::MatrixXd> combinedPZ,
                               Eigen::Ref<Eigen::MatrixXd> combinedP0,
                               const std::vector<Eigen::Matrix4d> &lorentz);


    std::tuple<std::vector<double>, std::vector<double>>
    JacobiCoordinatesMatrix(const Eigen::MatrixXd &particles,
                            const ClusterParams &params);

    std::tuple<std::vector<double>, std::vector<double>>
    JacobiCoordinates(const ParticleArray &particles,
                      const ClusterParams &params);


    std::tuple<double, double, double, double, double, double>
    fourBodyJacobi(const Coal::Particle &p1, const Coal::Particle &p2,
                   const Coal::Particle &n1, const Coal::Particle &n2);
}// namespace Coal
