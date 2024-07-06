//
// Created by mafu on 4/4/2024.
//
#pragma once
#include "ThreadPool.h"
#include "method.h"
#include <optional>

namespace Coal {

    class ClusterCalculator {
    public:
        static Eigen::MatrixXd
        Coalescence(const EventsMap &allEvents, const ClusterParams &params,
                    ResParamsMap &resolution,
                    std::optional<unsigned int> max_threads = std::nullopt);

    private:
        static std::mutex mtx;
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
            Eigen::Map<Eigen::MatrixXd> targetParticles;
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

        static Eigen::MatrixXd processEventstest(const EventsMap &allEvents,
                                                 const ClusterParams &params,
                                                 ResParamsMap &resolution,
                                                 unsigned int num_threads);
        static void processSubCellTest(const std::vector<Eigen::MatrixXd> &subCell,
                                       const ClusterParams &params,
                                       Eigen::MatrixXd &targetParticles, long &usedRows,
                                       int currentIndex);

        static Eigen::MatrixXd processEventsV2(const EventsMap &allEvents,
                                               const ClusterParams &params,
                                               unsigned int num_threads,
                                               ResParamsMap &resolution);

        static std::pair<int, int> getjobRange(int total_size, int num_threads,
                                               int thread_id);

        static Eigen::MatrixXd processEvents(const EventsMap &allEvents,
                                             const ClusterParams &params,
                                             ResParamsMap &resolution,
                                             unsigned int num_threads);

        static void processSubCell(const std::vector<Eigen::MatrixXd> &subCell,
                                   const ClusterParams &params,
                                   Eigen::MatrixXd &targetParticles,
                                   ResParamsMap &resolution, int currentIndex);

        static Eigen::MatrixXd mainLoop(const std::vector<Eigen::MatrixXd> &MArray,
                                        const ClusterParams &params);

        static Eigen::MatrixXd mainLoopV1(const std::vector<Eigen::MatrixXd> &MArray,
                                          const ClusterParams &params,
                                          unsigned int num_threads);

        static void mainLoopV2(const std::vector<Eigen::MatrixXd> &MArray,
                               const ClusterParams &params, int start_idx, int endidx,
                               Eigen::MatrixXd &threadOutputs, long &usedRows);

        static std::vector<int> jumpValidLoop(const std::vector<int> &multiIndex,
                                              const std::vector<int> &counts,
                                              int jumpIndex);

        static std::string createKeyFromMultiIndex(const std::vector<int> &multiIndex,
                                                   int index);

        static void CheckingPortfolioValidity(
                const Eigen::MatrixXd &ParticlesList, const ClusterParams &params,
                const std::vector<int> &counts, const std::vector<int> &multIndex,
                std::vector<int> &jumpMultiIndex,
                std::unordered_map<std::string, double> &distanceCache,
                std::vector<bool> &isValidCombination);

        static void setMatrix(MatrixMemPool &temMatrix, const Eigen::MatrixXd &particles,
                              const Eigen::MatrixXd &lastParticles, int pdg);

        static void vectorizationWithLastArray(
                const ClusterParams &params, MatrixMemPool &tempMatrixs,
                std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
                        &lorentzMatrixs);
    };
}// namespace Coal
