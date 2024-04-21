//
// Created by mafu on 1/6/2024.
//

#pragma once
#include "config.h"
#include "fileloader.h"
#include "flow.h"
#include <Eigen/Dense>
#include <utility>
namespace Coal {

    bool checkFileExits(const std::string &filename,
                        const std::vector<std::string> &labels,
                        const std::string &fileType);

    std::string constructFilename(const std::string &filename,
                                  const std::string &fileType, const std::string &label);

    std::vector<double> linspace(double start, double end, int num);

    EventsMap selectEvents(const EventsMap &eventMap, const ClusterParams &params,
                           ResParamsMap &resolution, int currentIndex,
                           std::vector<int> &eventIDList);

    std::vector<Eigen::MatrixXd> selectParticles(const EventsMap &eventMap,
                                                 const ClusterParams &params);

    //martix boost
    void lorenzBoostMatrix(Eigen::MatrixXd &particles, double beta_x, double beta_y,
                           double beta_z);

    void boostToComMatrix(Eigen::MatrixXd &particles);

    void calculateLorentz(const Eigen::VectorXd &betaX, const Eigen::VectorXd &betaY,
                          const Eigen::VectorXd &betaZ,
                          std::vector<Eigen::Matrix4d> &lorentzMatrixs);

    void applyLorentzBoost(Eigen::Ref<Eigen::MatrixXd> combinedX,
                           Eigen::Ref<Eigen::MatrixXd> combinedY,
                           Eigen::Ref<Eigen::MatrixXd> combinedZ,
                           Eigen::Ref<Eigen::MatrixXd> combinedT,
                           Eigen::Ref<Eigen::MatrixXd> combinedPX,
                           Eigen::Ref<Eigen::MatrixXd> combinedPY,
                           Eigen::Ref<Eigen::MatrixXd> combinedPZ,
                           Eigen::Ref<Eigen::MatrixXd> combinedP0,
                           const std::vector<Eigen::Matrix4d> &lorentz);

    std::tuple<std::vector<double>, std::vector<double>>
    JacobiCoordinatesMatrix_test(const Eigen::MatrixXd &particles,
                                 const ClusterParams &params);

    std::tuple<std::vector<double>, std::vector<double>>
    JacobiCoordinatesMatrix(const Eigen::MatrixXd &particles,
                            const ClusterParams &params);

    std::tuple<Eigen::VectorXd, Eigen::VectorXd>
    JacobiCoordinatesMatrix_Vec(const Eigen::MatrixXd &particles,
                              const ClusterParams &params);

}// namespace Coal
