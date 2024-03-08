//
// Created by mafu on 1/7/2024.
//
#pragma once
#include "method.h"

inline std::mutex mtx;
namespace Coal {

    struct ComputationMatrices {
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
        Eigen::Map<Eigen::MatrixXd> x_jacobi;
        Eigen::Map<Eigen::MatrixXd> y_jacobi;
        Eigen::Map<Eigen::MatrixXd> z_jacobi;
        Eigen::Map<Eigen::MatrixXd> px_jacobi;
        Eigen::Map<Eigen::MatrixXd> py_jacobi;
        Eigen::Map<Eigen::MatrixXd> pz_jacobi;
        Eigen::Map<Eigen::MatrixXd> dr;
        Eigen::Map<Eigen::MatrixXd> dp;
        Eigen::Map<Eigen::VectorXd> beta_x;
        Eigen::Map<Eigen::VectorXd> beta_y;
        Eigen::Map<Eigen::VectorXd> beta_z;
        Eigen::Map<Eigen::VectorXd> diff;

        // Constructor
        ComputationMatrices(double *memPool, int size_last, int N);
        void reset();
    };


    void clusterMatrix(const EventsMap &allEvents,
                       const std::string &outputFilename,
                       const std::string &ptFilename,
                       const ClusterParams &params, const YAML::Node &output);

    void clusterNoThread(const EventsMap &allEvents,
                         const std::string &outputFilename,
                         const std::string &ptFilename,
                         const ClusterParams &params,
                         const RapidityArray &rapidityArray);

    void clusterMix(const EventsMap &allEvents,
                    const std::string &outputFilename,
                    const std::string &ptFilename, const ClusterParams &params,
                    const RapidityArray &rapidityArray);

    void clusterThreadMatrix(const EventsMap &allEvents,
                             const std::string &outputFilename,
                             const std::string &ptFilename,
                             const ClusterParams &params,
                             const YAML::Node &output);
    void resetGlobalState();

    void clusterThreadV2(const EventsMap &allEvents,
                         const std::string &outputFilename,
                         const std::string &ptFilename,
                         const ClusterParams &params, const YAML::Node &config);

    void createThreadPool(unsigned int num_threads, const ClusterParams &params,
                          Eigen::MatrixXd &targetParticles);

    void processSubCell(const MultiParticleArray &subCell,
                        const ClusterParams &params,
                        ParticleArray &targetParticles);
    void processSubCellMatrix(const std::vector<Eigen::MatrixXd> &subCell,
                              const ClusterParams &params,
                              Eigen::MatrixXd &targetParticles);

    void outputTargetParticleMatrix(
            const Eigen::MatrixXd &targetParticles,
            std::ofstream &outputCluster,
            std::map<std::pair<double, double>, std::vector<double>> &ptArray,
            std::map<std::pair<double, double>, double> &yeildArray,
            const ClusterParams &params, const RapidityArray &rapidityArray,
            const bool &extended);

    void outputTargetParticle(
            const ParticleArray &targetParticles, std::ofstream &outputCluster,
            std::map<std::pair<double, double>, std::vector<double>> &ptArray,
            std::map<std::pair<double, double>, double> &yeildArray,
            const ClusterParams &params, const RapidityArray &rapidityArray);

    void outputTargetParticleBinary(
            const ParticleArray &targetParticles, std::ofstream &outputCluster,
            std::map<std::pair<double, double>, std::vector<double>> &ptArray,
            const ClusterParams &params, const RapidityArray &rapidityArray);

    std::vector<int> indexToMultiIndex(long long index,
                                       const std::vector<int> &counts);

    auto multiIndexToIndex(const std::vector<int> &multiIndex,
                           const std::vector<int> &counts) -> long long;

    bool incrementIndex(std::vector<int> &multiIndex,
                        const std::vector<int> &counts);

    ParticleArray mainLoopV2(const MultiParticleArray &ReactionParticles,
                             const ClusterParams &params);

    Eigen::MatrixXd mainLoopMatrix(const std::vector<Eigen::MatrixXd> &MArray,
                                   const ClusterParams &params);

    ParticleArray mainLoop(const MultiParticleArray &ReactionParticles,
                           const ClusterParams &params);

    //jump vaild loop
    std::vector<int> jumpValidLoop(const std::vector<int> &multiIndex,
                                   const std::vector<int> &counts,
                                   int jumpIndex);

    void nBodyCoal(const ParticleArray &Particles, Particle &targetParticle,
                   const ClusterParams &params, const std::vector<int> &counts,
                   const std::vector<int> &multIndex,
                   std::vector<int> &jumpMultiIndex,
                   std::vector<int> &lastMultiIndex);

    void nBodyCoalV2(const ParticleArray &Particles, Particle &targetParticle,
                     const ClusterParams &params,
                     const std::vector<int> &counts,
                     const std::vector<int> &multIndex,
                     std::vector<int> &jumpMultiIndex,
                     std::unordered_map<std::string, double> &distanceCache);

    std::string createKeyFromMultiIndex(const std::vector<int> &multiIndex,
                                        int index);

    std::vector<bool> checkCombinationList(
            const Eigen::MatrixXd &ParticlesList, const ClusterParams &params,
            const std::vector<int> &counts, const std::vector<int> &multIndex,
            std::vector<int> &jumpMultiIndex,
            std::unordered_map<std::string, double> &distanceCache);
    void batchProcessLastParticles(
            const Eigen::MatrixXd &Particles,
            const Eigen::MatrixXd &LastParticles, const ClusterParams &params,
            Eigen::Matrix<double, Eigen::Dynamic, 11> &targetParticles);

    void batchProcessLastParticlesCols(
            const Eigen::MatrixXd &Particles,
            const Eigen::MatrixXd &LastParticles, const ClusterParams &params,
            Eigen::Matrix<double, Eigen::Dynamic, 11> &targetParticles,
            ComputationMatrices &tempMatrices);

}// namespace Coal
