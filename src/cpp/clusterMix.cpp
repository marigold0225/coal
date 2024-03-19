//
// Created by mafu on 1/7/2024.
//
#include "../headers/clusterMix.h"
#include <condition_variable>
#include <execution>
#include <fstream>
#include <iostream>
#include <queue>
#include <ranges>
#include <thread>
#include <unordered_map>

std::vector<std::thread> workers;
std::mutex queue_mutex;
std::condition_variable cv;
std::queue<std::vector<Eigen::MatrixXd>> taskQueue;
bool stop_all = false;
std::mutex mtx;

Coal::ComputationMatrices::ComputationMatrices(double *memPool,
                                               const long size_last,
                                               const int N)
    : combinedX(memPool, size_last, N),
      combinedY(memPool += size_last * N, size_last, N),
      combinedZ(memPool += size_last * N, size_last, N),
      combinedT(memPool += size_last * N, size_last, N),
      combinedPX(memPool += size_last * N, size_last, N),
      combinedPY(memPool += size_last * N, size_last, N),
      combinedPZ(memPool += size_last * N, size_last, N),
      combinedP0(memPool += size_last * N, size_last, N),
      combinedMass(memPool += size_last * N, size_last, N),
      combinedProbability(memPool += size_last * N, size_last, N),
      x_jacobi(memPool += size_last * N, size_last, N - 1),
      y_jacobi(memPool += size_last * (N - 1), size_last, N - 1),
      z_jacobi(memPool += size_last * (N - 1), size_last, N - 1),
      px_jacobi(memPool += size_last * (N - 1), size_last, N - 1),
      py_jacobi(memPool += size_last * (N - 1), size_last, N - 1),
      pz_jacobi(memPool += size_last * (N - 1), size_last, N - 1),
      dr(memPool += size_last * (N - 1), size_last, N - 1),
      dp(memPool += size_last * (N - 1), size_last, N - 1),
      beta_x(memPool += size_last * (N - 1), size_last),
      beta_y(memPool += size_last, size_last),
      beta_z(memPool += size_last, size_last),
      diff(memPool += size_last, size_last) {}
void Coal::ComputationMatrices::reset() {
    combinedX.setZero();
    combinedY.setZero();
    combinedZ.setZero();
    combinedT.setZero();
    combinedPX.setZero();
    combinedPY.setZero();
    combinedPZ.setZero();
    combinedP0.setZero();
    combinedMass.setZero();
    combinedProbability.setZero();
    x_jacobi.setZero();
    y_jacobi.setZero();
    z_jacobi.setZero();
    px_jacobi.setZero();
    py_jacobi.setZero();
    pz_jacobi.setZero();
    dr.setZero();
    dp.setZero();
    beta_x.setZero();
    beta_y.setZero();
    beta_z.setZero();
    diff.setZero();
}

void Coal::clusterMatrix(const EventsMap &allEvents,
                         const std::string &outputFilename,
                         const std::string &ptFilename,
                         const ClusterParams &params,
                         const YAML::Node &output) {
    std::ofstream clusterOutput(outputFilename);
    std::ofstream ptOutput(ptFilename);
    std::map<RapidityRange, std::vector<double>> ptArray;
    std::map<RapidityRange, double> yeildArray;
    std::map<RapidityRange, std::vector<double>> v2Array;
    std::map<RapidityRange, std::vector<double>> countArray;
    const auto rapidityArray = output["RapidityRange"].as<RapidityArray>();
    const bool extended      = output["Extended"].as<bool>();
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(params.ptBins.second, 0.0);
        v2Array[rap]    = std::vector(params.ptBins.second, 0.0);
        countArray[rap] = std::vector(params.ptBins.second, 0.0);
        yeildArray[rap] = 0.0;
    }

    Eigen::MatrixXd targetParticles{};
    Eigen::MatrixXd resultParticles{};
    auto start_time = std::chrono::steady_clock::now();

    for (auto i = 0; i < params.Loop; ++i) {
        auto Cell    = selectEvents(allEvents, params);
        auto subCell = selectParticlesMatrix(Cell, params);
        if (subCell.empty()) {
            continue;
        }

        resultParticles = mainLoopMatrix(subCell, params);
        targetParticles.conservativeResize(
                targetParticles.rows() + resultParticles.rows(), 11);
        targetParticles.bottomRows(resultParticles.rows()) = resultParticles;
    }
    outputTargetParticleMatrix(targetParticles, clusterOutput, ptArray,
                               yeildArray, params, rapidityArray, extended,
                               v2Array, countArray);

    clusterOutput.close();

    for (const auto &rap: rapidityArray) {
        ptOutput << "Rapidity range: " << rap.first << "<y<" << rap.second
                 << ", cluster yield:"
                 << yeildArray[rap] / (rap.second - rap.first) << "\n";
        for (auto i = 0; i < params.ptBins.second; ++i) {
            double pt = params.ptBins.first / 2 +
                        static_cast<double>(i) * params.ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * params.ptBins.first *
                                std::abs((rap.second - rap.first)));

            v2Array[rap][i] /= countArray[rap][i];

            ptOutput << pt << " " << ptArray[rap][i] << " " << v2Array[rap][i]
                     << "\n";
        }
        ptOutput << "\n";
    }
    ptOutput.close();
    auto end_time = std::chrono::steady_clock::now();
    std::cout << "Time used: "
              << std::chrono::duration_cast<std::chrono::seconds>(end_time -
                                                                  start_time)
                         .count()

              << " s" << std::endl;
}

void Coal::resetGlobalState() {
    {
        std::unique_lock lock(queue_mutex);
        stop_all = true;
    }
    cv.notify_all();

    for (auto &worker: workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
    workers.clear();
    {
        std::lock_guard lock(queue_mutex);
        while (!taskQueue.empty()) {
            taskQueue.pop();
        }
    }
}


void Coal::clusterThreadV2(const EventsMap &allEvents,
                           const std::string &outputFilename,
                           const std::string &ptFilename,
                           const ClusterParams &params,
                           const YAML::Node &config) {
    resetGlobalState();

    {
        std::lock_guard lock(queue_mutex);
        stop_all = false;
    }

    std::ofstream clusterOutput(outputFilename);
    std::ofstream ptOutput(ptFilename);
    std::map<RapidityRange, std::vector<double>> ptArray;
    std::map<RapidityRange, double> yeildArray;
    std::map<RapidityRange, std::vector<double>> v2Array;
    std::map<RapidityRange, std::vector<double>> countArray;
    auto output              = config["Output"];
    const auto rapidityArray = output["RapidityRange"].as<RapidityArray>();
    const bool extended      = output["Extended"].as<bool>();
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(params.ptBins.second, 0.0);
        v2Array[rap]    = std::vector(params.ptBins.second, 0.0);
        countArray[rap] = std::vector(params.ptBins.second, 0.0);
        yeildArray[rap] = 0.0;
    }

    Eigen::MatrixXd targetParticles{};
    auto max_threads = config["General"]["Parallel"]["Cores"].as<unsigned>();
    std::cout << "max threads:" << max_threads << std::endl;
    if (max_threads > 0)
        max_threads -= 1;

    createThreadPool(max_threads, params, targetParticles);

    auto start_time = std::chrono::steady_clock::now();

    for (auto i = 0; i < params.Loop; ++i) {
        auto Cell    = selectEvents(allEvents, params);
        auto subCell = selectParticlesMatrix(Cell, params);
        if (subCell.empty()) {
            continue;
        }
        {
            std::lock_guard lock(queue_mutex);
            taskQueue.push(subCell);
        }
        cv.notify_one();
    }
    {
        std::lock_guard lock(queue_mutex);
        stop_all = true;
    }
    cv.notify_all();

    for (auto &thread: workers) {
        if (thread.joinable()) {
            thread.join();
        }
    }

    outputTargetParticleMatrix(targetParticles, clusterOutput, ptArray,
                               yeildArray, params, rapidityArray, extended,
                               v2Array, countArray);

    clusterOutput.close();

    for (const auto &rap: rapidityArray) {
        ptOutput << "Rapidity range: " << rap.first << "<y<" << rap.second
                 << ", cluster yield:"
                 << yeildArray[rap] / (rap.second - rap.first) << "\n";
        for (auto i = 0; i < params.ptBins.second; ++i) {
            double pt = params.ptBins.first / 2 +
                        static_cast<double>(i) * params.ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * params.ptBins.first *
                                std::abs((rap.second - rap.first)));

            v2Array[rap][i] = countArray[rap][i] > 0
                                      ? v2Array[rap][i] / countArray[rap][i]
                                      : 0.0;

            ptOutput << pt << " " << ptArray[rap][i] << " " << v2Array[rap][i]
                     << "\n";
        }
        ptOutput << "\n";
    }

    ptOutput.close();
    auto end_time = std::chrono::steady_clock::now();
    std::cout << "Time used: "
              << std::chrono::duration_cast<std::chrono::seconds>(end_time -
                                                                  start_time)
                         .count()

              << " s" << std::endl;
}
void Coal::createThreadPool(const unsigned num_threads,
                            const ClusterParams &params,
                            Eigen::MatrixXd &targetParticles) {
    for (unsigned i = 0; i < num_threads; ++i) {
        workers.emplace_back([&]() {
            while (true) {
                std::vector<Eigen::MatrixXd> subCell;
                {
                    std::unique_lock lock(queue_mutex);
                    cv.wait(lock,
                            [&]() { return stop_all || !taskQueue.empty(); });
                    if (stop_all && taskQueue.empty()) {
                        return;
                    }
                    subCell = std::move(taskQueue.front());
                    taskQueue.pop();
                }
                processSubCellMatrix(subCell, params, targetParticles);
            }
        });
    }
}

void Coal::processSubCellMatrix(const std::vector<Eigen::MatrixXd> &subCell,
                                const ClusterParams &params,
                                Eigen::MatrixXd &targetParticles) {
    const auto resultParticles = mainLoopMatrix(subCell, params);
    std::lock_guard lock(mtx);
    targetParticles.conservativeResize(
            targetParticles.rows() + resultParticles.rows(), 11);
    targetParticles.bottomRows(resultParticles.rows()) = resultParticles;
}


void Coal::outputTargetParticleMatrix(
        const Eigen::MatrixXd &targetParticles, std::ofstream &outputCluster,
        std::map<RapidityRange, std::vector<double>> &ptArray,
        std::map<RapidityRange, double> &yeildArray,
        const ClusterParams &params, const RapidityArray &rapidityArray,
        const bool &extended,
        std::map<RapidityRange, std::vector<double>> &v2Array,
        std::map<RapidityRange, std::vector<double>> &countArray) {
    const double d_pt        = params.ptBins.first;
    const int ptBins         = params.ptBins.second;
    const double totalEvents = params.eventFactor;
    outputCluster << "#events:" << totalEvents << " "
                  << "size:" << targetParticles.rows() << "\n";
    for (int i = 0; i < targetParticles.rows(); ++i) {

        const double pt = sqrt(targetParticles(i, 1) * targetParticles(i, 1) +
                               targetParticles(i, 2) * targetParticles(i, 2));

        const double v2 = (targetParticles(i, 1) * targetParticles(i, 1) -
                           targetParticles(i, 2) * targetParticles(i, 2)) /
                          (targetParticles(i, 1) * targetParticles(i, 1) +
                           targetParticles(i, 2) * targetParticles(i, 2));

        const double rapidity =
                0.5 * log((targetParticles(i, 4) + targetParticles(i, 3)) /
                          (targetParticles(i, 4) - targetParticles(i, 3)));
        for (const auto &rap: rapidityArray) {
            if (rapidity > rap.first && rapidity <= rap.second) {
                if (const int index = static_cast<int>(pt / d_pt);
                    index < ptBins && index >= 0) {
                    ptArray[rap][index] += targetParticles(i, 10) / totalEvents;
                    v2Array[rap][index] += v2 * targetParticles(i, 10);
                    countArray[rap][index] += targetParticles(i, 10);
                }
                yeildArray[rap] += targetParticles(i, 10) / totalEvents;
            }
        }
    }
    if (extended) {
        outputCluster << targetParticles << "\n";
    }
}

void Coal::outputTargetParticleBinary(
        const ParticleArray &targetParticles, std::ofstream &outputCluster,
        std::map<std::pair<double, double>, std::vector<double>> &ptArray,
        const ClusterParams &params, const RapidityArray &rapidityArray) {
    const double d_pt      = params.ptBins.first;
    const int ptBins       = params.ptBins.second;
    const int totalEvents  = static_cast<int>(params.eventFactor);
    const int particleSize = static_cast<int>(targetParticles.size());
    std::cout << "totalEvents:" << totalEvents << "\n";
    std::cout << "particleSize:" << particleSize << std::endl;
    outputCluster.write(reinterpret_cast<const char *>(&totalEvents),
                        sizeof(totalEvents));
    outputCluster.write(reinterpret_cast<const char *>(&particleSize),
                        sizeof(particleSize));
    for (auto &cluster: targetParticles) {
        const double pt       = cluster.getPT();
        const double rapidity = cluster.getRapidity();
        for (const auto &rap: rapidityArray) {
            if (rapidity > rap.first && rapidity < rap.second) {
                if (const int index = static_cast<int>(pt / d_pt);
                    index < ptBins) {
                    ptArray[rap][index] += cluster.probability / totalEvents;
                }
            }
        }
        outputCluster.write(reinterpret_cast<const char *>(&cluster),
                            sizeof(cluster));
    }
    outputCluster.close();
}

std::vector<int> Coal::indexToMultiIndex(const long long index,
                                         const std::vector<int> &counts) {
    std::vector multiIndex(counts.size(), 0);
    long long product = 1;
    std::ranges::for_each(counts, [&product](const int count) {
        if (count != 0) {
            product *= count;
        }
    });

    long long remaining = index;
    for (int i = 0; i < counts.size(); ++i) {
        if (counts[i] != 0) {
            product /= counts[i];
            multiIndex[i] = static_cast<int>(remaining / product);
            remaining     = remaining % product;
        } else {
            multiIndex[i] = 0;
        }
    }
    return multiIndex;
}

long long Coal::multiIndexToIndex(const std::vector<int> &multiIndex,
                                  const std::vector<int> &counts) {
    long long index   = 0;
    long long product = 1;
    std::ranges::for_each(counts, [&product](const int count) {
        if (count != 0) {
            product *= count;
        }
    });

    for (int i = 0; i < counts.size(); ++i) {
        if (counts[i] != 0) {
            product /= counts[i];
            index += multiIndex[i] * product;
        }
    }
    return index;
}

bool Coal::incrementIndex(std::vector<int> &multiIndex,
                          const std::vector<int> &counts) {
    const int dimension = static_cast<int>(counts.size());
    for (int i = dimension - 1; i >= 0; --i) {
        if (++multiIndex[i] < counts[i]) {
            return true;
        }
        multiIndex[i] = 0;
    }
    return false;
}

Eigen::MatrixXd Coal::mainLoopMatrix(const std::vector<Eigen::MatrixXd> &MArray,
                                     const ClusterParams &params) {
    const auto N         = params.NBody;
    const double factor  = params.eventFactor;
    const auto size_last = MArray[N - 1].rows();
    const long MAX_SIZE  = (size_last * N * 10)//combined Matrix
                          + (size_last * 4)    //beta_x, beta_y, beta_z, diff
                          + (size_last * (N - 1) * 8);
    std::vector memPool(MAX_SIZE, 0.0);
    const auto memPtr = memPool.data();

    ComputationMatrices tempMatrixs(memPtr, size_last, N);

    std::vector counts(N, 0);
    std::vector maxIndex(N, 0);
    for (auto i = 0; i < N; ++i) {
        counts[i]   = static_cast<int>(MArray[i].rows());
        maxIndex[i] = counts[i] - 1;
    }
    Eigen::MatrixXd targetParticleArray;
    std::unordered_map<std::string, double> distance_cache;
    double yieldAll    = 0.0;
    double yieldSelect = 0.0;
    long long loop     = 0;
    std::vector multiIndex(N, 0);
    Eigen::MatrixXd particles(N, 11);
    Eigen::Matrix<double, Eigen::Dynamic, 11> targetParticle(size_last, 11);
    while (true) {
        auto jumpMultiIndex = multiIndex;
        for (auto j = 0; j < N; ++j) {
            particles.row(j) = MArray[j].row(multiIndex[j]);
        }

        auto isValidCombination =
                checkCombinationList(particles, params, counts, multiIndex,
                                     jumpMultiIndex, distance_cache);
        if (std::ranges::all_of(isValidCombination,
                                [](const bool x) { return x; })) {

            targetParticle.setZero();
            tempMatrixs.reset();
            batchProcessLastParticlesCols(
                    particles.block(0, 0, particles.rows() - 1,
                                    particles.cols()),
                    MArray[N - 1], params, targetParticle, tempMatrixs);
            yieldAll += targetParticle.col(10).sum();
            if (Eigen::Array<bool, Eigen::Dynamic, 1> mask =
                        targetParticle.col(10).array() > params.probabilityCut;
                mask.any()) {
                std::vector<int> indices;
                for (int i = 0; i < mask.size(); ++i) {
                    if (mask(i)) {
                        indices.push_back(i);
                    }
                }

                Eigen::MatrixXd filteredTargetParticle(indices.size(),
                                                       targetParticle.cols());
                for (long i = 0; i < indices.size(); ++i) {
                    filteredTargetParticle.row(i) =
                            targetParticle.row(indices[i]);
                }

                yieldSelect += filteredTargetParticle.col(10).sum();

                targetParticleArray.conservativeResize(
                        targetParticleArray.rows() +
                                filteredTargetParticle.rows(),
                        11);
                targetParticleArray.bottomRows(filteredTargetParticle.rows()) =
                        filteredTargetParticle;
            }

            multiIndex = jumpValidLoop(multiIndex, counts, N - 3);
            if (multiIndex == maxIndex) {
                break;
            }

        } else {
            if (jumpMultiIndex == maxIndex) {
                break;
            }
            multiIndex = jumpMultiIndex;
            continue;
        }
        loop++;
    }

    targetParticleArray.col(10) *= yieldAll / yieldSelect;
    std::cout << "loop:" << loop << " "
              << "yieldAll:" << yieldAll << " "
              << "yieldSelect:" << yieldSelect << " "
              << "factor:" << factor << " "
              << "size:" << targetParticleArray.rows() << "\n";
    return targetParticleArray;
}


// Coal::ParticleArray Coal::mainLoop(const MultiParticleArray &ReactionParticles,
//                                    const ClusterParams &params) {
//     const auto N        = params.NBody;
//     const double factor = params.eventFactor;
//     std::vector counts(N, 0);
//     long long loopMax = 1;
//     for (auto i = 0; i < N; ++i) {
//         counts[i] = static_cast<int>(ReactionParticles[i].size());
//         loopMax *= counts[i];
//     }
//     ParticleArray targetParticleArray{};
//     double yieldAll    = 0.0;
//     double yieldSelect = 0.0;
//     long long loop     = 0;
//     std::vector lastmultiIndex(N, -1);
//     for (long long i = 0; i < loopMax;) {
//
//         const auto multiIndex = indexToMultiIndex(i, counts);
//         auto jumpMultiIndex   = multiIndex;
//         ParticleArray particles{};
//         for (auto j = 0; j < N; ++j) {
//             particles.push_back(ReactionParticles[j][multiIndex[j]]);
//         }
//         Particle targetParticle{};
//
//         nBodyCoal(particles, targetParticle, params, counts, multiIndex,
//                   jumpMultiIndex, lastmultiIndex);
//
//         if (jumpMultiIndex != multiIndex) {
//             i = multiIndexToIndex(jumpMultiIndex, counts);
//             if (i >= loopMax) {
//                 break;
//             }
//             continue;
//         }
//         yieldAll += targetParticle.probability;
//
//         if (targetParticle.probability > 0.0 &&
//             targetParticle.probability > params.probabilityCut) {
//             yieldSelect += targetParticle.probability;
//             targetParticleArray.push_back(targetParticle);
//         }
//         loop++;
//         ++i;
//     }
//     for (auto &particle: targetParticleArray) {
//         particle.probability *= yieldAll / yieldSelect;
//     }
//     std::cout << "loopMax:" << loopMax << " "
//               << "loopCut:" << loop << " "
//               << "yieldAll:" << yieldAll << " "
//               << "yieldSelect:" << yieldSelect << " "
//               << "factor:" << factor << " "
//               << "size:" << targetParticleArray.size() << "\n";
//
//     return targetParticleArray;
// }

std::vector<int> Coal::jumpValidLoop(const std::vector<int> &multiIndex,
                                     const std::vector<int> &counts,
                                     const int jumpIndex) {
    if (jumpIndex >= static_cast<int>(multiIndex.size()) - 1) {
        return multiIndex;
    }
    auto newMultiIndex = multiIndex;
    newMultiIndex[jumpIndex + 1]++;

    for (int i = jumpIndex + 1; i >= 0; --i) {
        if (newMultiIndex[i] >= counts[i]) {
            newMultiIndex[i] = 0;
            if (i > 0) {
                newMultiIndex[i - 1]++;
            }
        } else {
            break;
        }
    }
    for (int i = jumpIndex + 2; i < newMultiIndex.size(); ++i) {
        newMultiIndex[i] = 0;
    }

    if (std::ranges::all_of(newMultiIndex,
                            [](const int x) { return x == 0; })) {
        for (int i = 0; i < newMultiIndex.size(); ++i) {
            newMultiIndex[i] = counts[i] - 1;
        }
    }

    return newMultiIndex;
}

std::string Coal::createKeyFromMultiIndex(const std::vector<int> &multiIndex,
                                          const int index) {
    std::string key;
    for (auto i = 0; i <= index; ++i) {
        key += std::to_string(multiIndex[i]);
        if (i < index) {
            key += "-";
        }
    }
    return key;
}
std::vector<bool> Coal::checkCombinationList(
        const Eigen::MatrixXd &ParticlesList, const ClusterParams &params,
        const std::vector<int> &counts, const std::vector<int> &multIndex,
        std::vector<int> &jumpMultiIndex,
        std::unordered_map<std::string, double> &distanceCache) {
    constexpr double hbar  = 0.19733;
    constexpr double hbar2 = hbar * hbar;
    const auto N           = static_cast<int>(ParticlesList.rows());
    if (N < 3) {
        return std::vector(N, true);
    }
    std::vector isValidCombination(N - 2, true);

    for (auto i = 1; i < N - 1; ++i) {
        auto dis_temp   = 0.0;
        std::string key = createKeyFromMultiIndex(multIndex, i);

        if (auto it = distanceCache.find(key); it != distanceCache.end()) {
            dis_temp = it->second;
        } else {

            for (auto index = distanceCache.begin();
                 index != distanceCache.end();) {
                if (std::ranges::count(index->first, '-') == i) {
                    index = distanceCache.erase(index);
                } else {
                    ++index;
                }
            }

            Eigen::MatrixXd tempParticle = ParticlesList.topRows(i + 1);
            const auto boost_particles   = boostToComMatrix(tempParticle);

            auto [diff_r, diff_p] =
                    JacobiCoordinatesMatrix(boost_particles, params);
            for (auto j = 0; j < i; ++j) {
                dis_temp += (diff_r[j] * diff_r[j] / params.SigArray[j] /
                                     params.SigArray[j] +
                             diff_p[j] * diff_p[j] * params.SigArray[j] *
                                     params.SigArray[j] / hbar2);
            }
            distanceCache[key] = dis_temp;
        }

        if (dis_temp > (params.precision + (i - 1) * 8)) {
            isValidCombination[i - 1] = false;
            jumpMultiIndex            = jumpValidLoop(multIndex, counts, i - 1);
            return isValidCombination;
        }
    }
    return isValidCombination;
}
void Coal::batchProcessLastParticlesCols(
        const Eigen::MatrixXd &Particles, const Eigen::MatrixXd &LastParticles,
        const ClusterParams &params,
        Eigen::Matrix<double, -1, 11> &targetParticles,
        ComputationMatrices &tempMatrices) {
    constexpr double hbar  = 0.19733;
    constexpr double hbar2 = hbar * hbar;
    const auto size        = static_cast<int>(Particles.rows());
    const auto size_last   = static_cast<int>(LastParticles.rows());
    const auto N           = params.NBody;

    for (auto i = 0; i < size; ++i) {
        tempMatrices.combinedX.col(i).setConstant(Particles(i, 6));

        tempMatrices.combinedY.col(i).setConstant(Particles(i, 7));

        tempMatrices.combinedZ.col(i).setConstant(Particles(i, 8));
        tempMatrices.combinedT.col(i).setConstant(Particles(i, 9));

        tempMatrices.combinedPX.col(i).setConstant(Particles(i, 1));

        tempMatrices.combinedPY.col(i).setConstant(Particles(i, 2));

        tempMatrices.combinedPZ.col(i).setConstant(Particles(i, 3));

        tempMatrices.combinedP0.col(i).setConstant(Particles(i, 5));

        tempMatrices.combinedMass.col(i).setConstant(Particles(i, 4));

        tempMatrices.combinedProbability.col(i).setConstant(Particles(i, 10));
    }
    tempMatrices.combinedX.col(N - 1)           = LastParticles.col(6);
    tempMatrices.combinedY.col(N - 1)           = LastParticles.col(7);
    tempMatrices.combinedZ.col(N - 1)           = LastParticles.col(8);
    tempMatrices.combinedT.col(N - 1)           = LastParticles.col(9);
    tempMatrices.combinedPX.col(N - 1)          = LastParticles.col(1);
    tempMatrices.combinedPY.col(N - 1)          = LastParticles.col(2);
    tempMatrices.combinedPZ.col(N - 1)          = LastParticles.col(3);
    tempMatrices.combinedP0.col(N - 1)          = LastParticles.col(5);
    tempMatrices.combinedMass.col(N - 1)        = LastParticles.col(4);
    tempMatrices.combinedProbability.col(N - 1) = LastParticles.col(10);

    targetParticles.col(0) = params.pdg * Eigen::VectorXd::Ones(size_last);
    targetParticles.col(1) = tempMatrices.combinedPX.rowwise().sum();
    targetParticles.col(2) = tempMatrices.combinedPY.rowwise().sum();
    targetParticles.col(3) = tempMatrices.combinedPZ.rowwise().sum();
    targetParticles.col(4) = tempMatrices.combinedP0.rowwise().sum();
    targetParticles.col(5) = tempMatrices.combinedMass.rowwise().sum();
    targetParticles.col(6) = tempMatrices.combinedX.rowwise().mean();
    targetParticles.col(7) = tempMatrices.combinedY.rowwise().mean();
    targetParticles.col(8) = tempMatrices.combinedZ.rowwise().mean();
    targetParticles.col(9) = tempMatrices.combinedT.rowwise().maxCoeff();

    tempMatrices.beta_x =
            targetParticles.col(1).array() / targetParticles.col(4).array();
    tempMatrices.beta_y =
            targetParticles.col(2).array() / targetParticles.col(4).array();
    tempMatrices.beta_z =
            targetParticles.col(3).array() / targetParticles.col(4).array();

    const auto lambda = calculateLorentz(
            tempMatrices.beta_x, tempMatrices.beta_y, tempMatrices.beta_z);
    applyLorentzBoost(tempMatrices.combinedX, tempMatrices.combinedY,
                      tempMatrices.combinedZ, tempMatrices.combinedT,
                      tempMatrices.combinedPX, tempMatrices.combinedPY,
                      tempMatrices.combinedPZ, tempMatrices.combinedP0, lambda);

    tempMatrices.combinedX =
            (tempMatrices.combinedT.rowwise()
                     .maxCoeff()
                     .replicate(1, tempMatrices.combinedT.cols())
                     .array() -
             tempMatrices.combinedT.array()) *
                    tempMatrices.combinedPX.array() /
                    tempMatrices.combinedP0.array() +
            tempMatrices.combinedX.array();
    tempMatrices.combinedY =
            (tempMatrices.combinedT.rowwise()
                     .maxCoeff()
                     .replicate(1, tempMatrices.combinedT.cols())
                     .array() -
             tempMatrices.combinedT.array()) *
                    tempMatrices.combinedPY.array() /
                    tempMatrices.combinedP0.array() +
            tempMatrices.combinedY.array();
    tempMatrices.combinedZ =
            (tempMatrices.combinedT.rowwise()
                     .maxCoeff()
                     .replicate(1, tempMatrices.combinedT.cols())
                     .array() -
             tempMatrices.combinedT.array()) *
                    tempMatrices.combinedPZ.array() /
                    tempMatrices.combinedP0.array() +
            tempMatrices.combinedZ.array();

    tempMatrices.x_jacobi =
            (params.M[N - 2] * tempMatrices.combinedX.transpose())
                    .transpose()
                    .rightCols(N - 1);
    tempMatrices.y_jacobi =
            (params.M[N - 2] * tempMatrices.combinedY.transpose())
                    .transpose()
                    .rightCols(N - 1);
    tempMatrices.z_jacobi =
            (params.M[N - 2] * tempMatrices.combinedZ.transpose())
                    .transpose()
                    .rightCols(N - 1);
    tempMatrices.px_jacobi =
            (params.M_inv_t[N - 2] * tempMatrices.combinedPX.transpose())
                    .transpose()
                    .rightCols(N - 1);
    tempMatrices.py_jacobi =
            (params.M_inv_t[N - 2] * tempMatrices.combinedPY.transpose())
                    .transpose()
                    .rightCols(N - 1);
    tempMatrices.pz_jacobi =
            (params.M_inv_t[N - 2] * tempMatrices.combinedPZ.transpose())
                    .transpose()
                    .rightCols(N - 1);


    tempMatrices.dr = tempMatrices.x_jacobi.array().square() +
                      tempMatrices.y_jacobi.array().square() +
                      tempMatrices.z_jacobi.array().square();
    tempMatrices.dp = tempMatrices.px_jacobi.array().square() +
                      tempMatrices.py_jacobi.array().square() +
                      tempMatrices.pz_jacobi.array().square();

    for (auto i = 0; i < N - 1; ++i) {
        tempMatrices.diff +=
                (tempMatrices.dr.col(i).array() / params.SigArray[i] /
                         params.SigArray[i] +
                 tempMatrices.dp.col(i).array() * params.SigArray[i] *
                         params.SigArray[i] / hbar2)
                        .matrix();
    }
    targetParticles.col(10) =
            params.gc * params.probFactor *
            tempMatrices.combinedProbability.rowwise().prod().array() *
            (-tempMatrices.diff.array()).exp();
}
