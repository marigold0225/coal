//
// Created by mafu on 1/7/2024.
//
#pragma once
#include "method.h"
namespace Coal {

    struct MatrixMemPool {

        Eigen::Map<Eigen::MatrixXd> combinedX;
        Eigen::Map<Eigen::MatrixXd> combinedY;
        Eigen::Map<Eigen::MatrixXd> combinedZ;
        Eigen::Map<Eigen::MatrixXd> combinedT;
        Eigen::Map<Eigen::MatrixXd> combinedPX;
        Eigen::Map<Eigen::MatrixXd> combinedPY;
        Eigen::Map<Eigen::MatrixXd> combinedPZ;
        Eigen::Map<Eigen::MatrixXd> combinedP0;
        Eigen::Map<Eigen::MatrixXd> combinedMass;
        Eigen::Map<Eigen::MatrixXd> combinedProbability;
        Eigen::Map<Eigen::MatrixXd> temReplicatedMaxCoeff;
        Eigen::Map<Eigen::MatrixXd> dr;
        Eigen::Map<Eigen::MatrixXd> dp;
        Eigen::Map<Eigen::VectorXd> beta_x;
        Eigen::Map<Eigen::VectorXd> beta_y;
        Eigen::Map<Eigen::VectorXd> beta_z;
        Eigen::Map<Eigen::VectorXd> diff;

        // Constructor
        MatrixMemPool(double *memPool, int size_last, int N);
        void reset();
    };


    void singleThreadedParallelism(const EventsMap &allEvents,
                                   const std::string &outputFilename,
                                   const std::string &ptFilename,
                                   const ClusterParams &params,
                                   const YAML::Node &output);

    void resetGlobalState();

    void multithreadedParallelism(const EventsMap &allEvents,
                                  const std::string &outputFilename,
                                  const std::string &ptFilename,
                                  const ClusterParams &params,
                                  const YAML::Node &config);

    void createThreadPool(unsigned int num_threads, const ClusterParams &params,
                          Eigen::MatrixXd &targetParticles);

    void processSubCell(const std::vector<Eigen::MatrixXd> &subCell,
                        const ClusterParams &params,
                        Eigen::MatrixXd &targetParticles);

    void outputMatrix(const Eigen::MatrixXd &targetParticles,
                      std::ofstream &outputCluster,
                      std::map<RapidityRange, std::vector<double>> &ptArray,
                      std::map<RapidityRange, double> &yeildArray,
                      const ClusterParams &params,
                      const RapidityArray &rapidityArray, const bool &extended,
                      std::map<RapidityRange, std::vector<double>> &v2Array,
                      std::map<RapidityRange, std::vector<double>> &countArray);

    // void outputV2(const Eigen::MatrixXd &targetParticles,
    //               std::map<RapidityRange, std::vector<double>> &v2Array,
    //               const ClusterParams &params,
    //               const RapidityArray &rapidityArray);

    void outputBinary(
            const ParticleArray &targetParticles, std::ofstream &outputCluster,
            std::map<std::pair<double, double>, std::vector<double>> &ptArray,
            const ClusterParams &params, const RapidityArray &rapidityArray);

    std::vector<int> indexToMultiIndex(long long index,
                                       const std::vector<int> &counts);

    auto multiIndexToIndex(const std::vector<int> &multiIndex,
                           const std::vector<int> &counts) -> long long;

    bool incrementIndex(std::vector<int> &multiIndex,
                        const std::vector<int> &counts);


    Eigen::MatrixXd mainLoop(const std::vector<Eigen::MatrixXd> &MArray,
                             const ClusterParams &params);


    //jump vaild loop
    std::vector<int> jumpValidLoop(const std::vector<int> &multiIndex,
                                   const std::vector<int> &counts,
                                   int jumpIndex);

    std::string createKeyFromMultiIndex(const std::vector<int> &multiIndex,
                                        int index);

    std::vector<bool> CheckingPortfolioValidity(
            const Eigen::MatrixXd &ParticlesList, const ClusterParams &params,
            const std::vector<int> &counts, const std::vector<int> &multIndex,
            std::vector<int> &jumpMultiIndex,
            std::unordered_map<std::string, double> &distanceCache);

    void setMatrix(MatrixMemPool &temMatrix, const Eigen::MatrixXd &particles,
                   const Eigen::MatrixXd &lastParticles);

    void vectorizationWithLastArray(
            int size_last, const ClusterParams &params,
            Eigen::Matrix<double, Eigen::Dynamic, 11> &targetParticles,
            MatrixMemPool &tempMatrixs,
            std::vector<Eigen::Matrix4d> &lorentzMatrixs);

}// namespace Coal
