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

Coal::MatrixMemPool::MatrixMemPool(double *memPool, const int size_last,
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
      temReplicatedMaxCoeff(memPool += size_last * N, size_last, N),
      dr(memPool += size_last * N, size_last, N - 1),
      dp(memPool += size_last * (N - 1), size_last, N - 1),
      beta_x(memPool += size_last * (N - 1), size_last),
      beta_y(memPool += size_last, size_last),
      beta_z(memPool += size_last, size_last),
      diff(memPool += size_last, size_last) {}

void Coal::MatrixMemPool::reset() {
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
    temReplicatedMaxCoeff.setZero();
    dr.setZero();
    dp.setZero();
    beta_x.setZero();
    beta_y.setZero();
    beta_z.setZero();
    diff.setZero();
}

void Coal::singleThreadedParallelism(const EventsMap &allEvents,
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

        resultParticles = mainLoop(subCell, params);
        targetParticles.conservativeResize(
                targetParticles.rows() + resultParticles.rows(), 11);
        targetParticles.bottomRows(resultParticles.rows()) = resultParticles;
    }
    outputMatrix(targetParticles, clusterOutput, ptArray, yeildArray, params,
                 rapidityArray, extended, v2Array, countArray);

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

            v2Array[rap][i] /= countArray[rap][i] > 0 ? countArray[rap][i] : 1;

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


void Coal::multithreadedParallelism(const EventsMap &allEvents,
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

    outputMatrix(targetParticles, clusterOutput, ptArray, yeildArray, params,
                 rapidityArray, extended, v2Array, countArray);

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
                processSubCell(subCell, params, targetParticles);
            }
        });
    }
}

void Coal::processSubCell(const std::vector<Eigen::MatrixXd> &subCell,
                          const ClusterParams &params,
                          Eigen::MatrixXd &targetParticles) {
    const auto resultParticles = mainLoop(subCell, params);
    std::lock_guard lock(mtx);
    targetParticles.conservativeResize(
            targetParticles.rows() + resultParticles.rows(), 11);
    targetParticles.bottomRows(resultParticles.rows()) = resultParticles;
}


void Coal::outputMatrix(
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

        const double rapidity =
                0.5 * log((targetParticles(i, 4) + targetParticles(i, 3)) /
                          (targetParticles(i, 4) - targetParticles(i, 3)));
        for (const auto &rap: rapidityArray) {
            if (rapidity > rap.first && rapidity <= rap.second) {
                if (const int index = static_cast<int>(pt / d_pt);
                    index < ptBins && index >= 0) {
                    ptArray[rap][index] += targetParticles(i, 10) / totalEvents;
                }
                yeildArray[rap] += targetParticles(i, 10) / totalEvents;
            }
        }
    }
    if (extended) {
        outputCluster << targetParticles << "\n";
    }
}
// void Coal::outputV2(const Eigen::MatrixXd &targetParticles,
//                     std::map<RapidityRange, std::vector<double>> &v2Array,
//                     const ClusterParams &params,
//                     const RapidityArray &rapidityArray) {
//     const double d_pt        = params.ptBins.first;
//     const int ptBins         = params.ptBins.second;
//     const double totalEvents = params.eventFactor;
//     double sin2phi           = 0.0;
//     double cos2phi           = 0.0;
//     double sin2phi_A         = 0.0;
//     double cos2phi_A         = 0.0;
//     double sin2phi_B         = 0.0;
//     double cos2phi_B         = 0.0;
//     std::vector cos_i(targetParticles.rows(), 0.0);
//     std::vector sin_i(targetParticles.rows(), 0.0);
//     std::vector pt_i(targetParticles.rows(), 0.0);
//     std::vector phi_i(targetParticles.rows(), 0.0);
//
//
//     for (int i = 0; i < targetParticles.rows(); ++i) {
//         const double pt = sqrt(targetParticles(i, 1) * targetParticles(i, 1) +
//                                targetParticles(i, 2) * targetParticles(i, 2));
//         const double phi =
//                 std::atan2(targetParticles(i, 2), targetParticles(i, 1));
//         const double p_total =
//                 sqrt(targetParticles(i, 1) * targetParticles(i, 1) +
//                      targetParticles(i, 2) * targetParticles(i, 2) +
//                      targetParticles(i, 3) * targetParticles(i, 3));
//         const double pesudoRapidity =
//                 0.5 * log((p_total + targetParticles(i, 3)) /
//                           (p_total - targetParticles(i, 3)));
//         cos_i[i] = cos(2 * phi);
//         sin_i[i] = sin(2 * phi);
//         pt_i[i]  = pt;
//         phi_i[i] = phi;
//
//         if (pesudoRapidity > 1.0 || pesudoRapidity < -1.0) {
//             continue;
//         }
//
//         sin2phi += targetParticles(i, 10) * sin(2 * phi) * pt;
//         cos2phi += targetParticles(i, 10) * cos(2 * phi) * pt;
//
//         if (pesudoRapidity > -1.0 && pesudoRapidity < -0.05) {
//             sin2phi_A += targetParticles(i, 10) * sin(2 * phi) * pt;
//             cos2phi_A += targetParticles(i, 10) * cos(2 * phi) * pt;
//         }
//         if (pesudoRapidity > 0.05 && pesudoRapidity < 1.0) {
//             sin2phi_B += targetParticles(i, 10) * sin(2 * phi) * pt;
//             cos2phi_B += targetParticles(i, 10) * cos(2 * phi) * pt;
//         }
//     }
//     double Phi_A    = std::atan2(sin2phi_A, cos2phi_A) / 2;
//     double Phi_B    = std::atan2(sin2phi_B, cos2phi_B) / 2;
//     double Cos_2phi = cos(2 * (Phi_A - Phi_B)) / totalEvents;
//     double Refsub   = std::sqrt(Cos_2phi);
//
//     double x       = solveEquation(Refsub);
//     double Reffull = getRef(x * sqrt(2));
//
//     for (int i = 0; i < targetParticles.rows(); ++i) {
//         double cos_w  = cos2phi - cos_i[i] * pt_i[i];
//         double sin_w  = sin2phi - sin_i[i] * pt_i[i];
//         double psi    = std::atan2(sin_w, cos_w) / 2;
//         double phi    = phi_i[i];
//         double v2_sub = cos(2 * (phi - psi));
//
//         double v2 = v2_sub / Reffull;
//
//         const auto rapidity =
//                 0.5 * log((targetParticles(i, 4) + targetParticles(i, 3)) /
//                           (targetParticles(i, 4) - targetParticles(i, 3)));
//         for (const auto &rap: rapidityArray) {
//             if (rapidity > rap.first && rapidity <= rap.second) {
//                 if (const int index = static_cast<int>(pt_i[i] / d_pt);
//                     index < ptBins && index >= 0) {
//                     v2Array[rap][index] += v2 * targetParticles(i, 10);
//                 }
//             }
//         }
//     }
// }

void Coal::outputBinary(
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

Eigen::MatrixXd Coal::mainLoop(const std::vector<Eigen::MatrixXd> &MArray,
                               const ClusterParams &params) {
    const auto N        = params.NBody;
    const double factor = params.eventFactor;
    const int size_last = static_cast<int>(MArray[N - 1].rows());

    const long MAX_SIZE = (size_last * N * 11)//combined Matrix
                          + (size_last * 4)   //beta_x, beta_y, beta_z, diff
                          + (size_last * (N - 1) * 2);
    std::vector memPool(MAX_SIZE, 0.0);
    const auto memPtr = memPool.data();

    MatrixMemPool tempMatrixs(memPtr, size_last, N);

    std::vector<Eigen::Matrix4d> lorentzMatrixs(size_last,
                                                Eigen::Matrix4d::Zero());
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
    // std::vector particleChanged(N, true);
    // std::vector lastParticleIndex(N, -1);
    while (true) {
        auto jumpMultiIndex = multiIndex;
        for (auto j = 0; j < N; ++j) {
            particles.row(j) = MArray[j].row(multiIndex[j]);
        }

        auto isValidCombination =
                CheckingPortfolioValidity(particles, params, counts, multiIndex,
                                          jumpMultiIndex, distance_cache);
        if (std::ranges::all_of(isValidCombination,
                                [](const bool x) { return x; })) {
            //
            // for (auto j = 0; j < N; ++j) {
            //     particleChanged[j] = lastParticleIndex[j] != multiIndex[j];
            // }
            targetParticle.setZero();
            tempMatrixs.reset();
            setMatrix(tempMatrixs,
                      particles.block(0, 0, particles.rows() - 1,
                                      particles.cols()),
                      MArray[N - 1]);
            vectorizationWithLastArray(size_last, params, targetParticle,
                                       tempMatrixs, lorentzMatrixs);
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
            // lastParticleIndex = multiIndex;
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
std::vector<bool> Coal::CheckingPortfolioValidity(
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
                    JacobiCoordinatesMatrix_test2(boost_particles, params);
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
void Coal::setMatrix(MatrixMemPool &temMatrix, const Eigen::MatrixXd &particles,
                     const Eigen::MatrixXd &lastParticles) {

    const int size = static_cast<int>(particles.rows());
    for (auto i = 0; i < size; ++i) {
        temMatrix.combinedX.col(i).setConstant(particles(i, 6));
        temMatrix.combinedY.col(i).setConstant(particles(i, 7));
        temMatrix.combinedZ.col(i).setConstant(particles(i, 8));
        temMatrix.combinedT.col(i).setConstant(particles(i, 9));
        temMatrix.combinedPX.col(i).setConstant(particles(i, 1));
        temMatrix.combinedPY.col(i).setConstant(particles(i, 2));
        temMatrix.combinedPZ.col(i).setConstant(particles(i, 3));
        temMatrix.combinedP0.col(i).setConstant(particles(i, 5));
        temMatrix.combinedMass.col(i).setConstant(particles(i, 4));
        temMatrix.combinedProbability.col(i).setConstant(particles(i, 10));
    }
    temMatrix.combinedX.col(size)           = lastParticles.col(6);
    temMatrix.combinedY.col(size)           = lastParticles.col(7);
    temMatrix.combinedZ.col(size)           = lastParticles.col(8);
    temMatrix.combinedT.col(size)           = lastParticles.col(9);
    temMatrix.combinedPX.col(size)          = lastParticles.col(1);
    temMatrix.combinedPY.col(size)          = lastParticles.col(2);
    temMatrix.combinedPZ.col(size)          = lastParticles.col(3);
    temMatrix.combinedP0.col(size)          = lastParticles.col(5);
    temMatrix.combinedMass.col(size)        = lastParticles.col(4);
    temMatrix.combinedProbability.col(size) = lastParticles.col(10);
}
void Coal::vectorizationWithLastArray(
        const int size_last, const ClusterParams &params,
        Eigen::Matrix<double, -1, 11> &targetParticles,
        MatrixMemPool &tempMatrixs,
        std::vector<Eigen::Matrix4d> &lorentzMatrixs) {

    constexpr double hbar  = 0.19733;
    constexpr double hbar2 = hbar * hbar;
    const auto N           = params.NBody;

    for (auto &matrix: lorentzMatrixs) {
        matrix.setZero();
    }

    targetParticles.col(0) = params.pdg * Eigen::VectorXd::Ones(size_last);
    targetParticles.col(1) = tempMatrixs.combinedPX.rowwise().sum();
    targetParticles.col(2) = tempMatrixs.combinedPY.rowwise().sum();
    targetParticles.col(3) = tempMatrixs.combinedPZ.rowwise().sum();
    targetParticles.col(4) = tempMatrixs.combinedP0.rowwise().sum();
    targetParticles.col(5) = tempMatrixs.combinedMass.rowwise().sum();
    targetParticles.col(6) = tempMatrixs.combinedX.rowwise().mean();
    targetParticles.col(7) = tempMatrixs.combinedY.rowwise().mean();
    targetParticles.col(8) = tempMatrixs.combinedZ.rowwise().mean();
    targetParticles.col(9) = tempMatrixs.combinedT.rowwise().maxCoeff();

    tempMatrixs.beta_x =
            targetParticles.col(1).array() / targetParticles.col(4).array();
    tempMatrixs.beta_y =
            targetParticles.col(2).array() / targetParticles.col(4).array();
    tempMatrixs.beta_z =
            targetParticles.col(3).array() / targetParticles.col(4).array();

    calculateLorentz(tempMatrixs.beta_x, tempMatrixs.beta_y, tempMatrixs.beta_z,
                     lorentzMatrixs);

    applyLorentzBoost(tempMatrixs.combinedX, tempMatrixs.combinedY,
                      tempMatrixs.combinedZ, tempMatrixs.combinedT,
                      tempMatrixs.combinedPX, tempMatrixs.combinedPY,
                      tempMatrixs.combinedPZ, tempMatrixs.combinedP0,
                      lorentzMatrixs);
    tempMatrixs.temReplicatedMaxCoeff =
            tempMatrixs.combinedT.rowwise()
                    .maxCoeff()
                    .replicate(1, tempMatrixs.combinedT.cols())
                    .array() -
            tempMatrixs.combinedT.array();

    tempMatrixs.combinedX = tempMatrixs.temReplicatedMaxCoeff.array() *
                                    tempMatrixs.combinedPX.array() /
                                    tempMatrixs.combinedP0.array() +
                            tempMatrixs.combinedX.array();
    tempMatrixs.combinedY = tempMatrixs.temReplicatedMaxCoeff.array() *
                                    tempMatrixs.combinedPY.array() /
                                    tempMatrixs.combinedP0.array() +
                            tempMatrixs.combinedY.array();
    tempMatrixs.combinedZ = tempMatrixs.temReplicatedMaxCoeff.array() *
                                    tempMatrixs.combinedPZ.array() /
                                    tempMatrixs.combinedP0.array() +
                            tempMatrixs.combinedZ.array();

    tempMatrixs.dr = ((params.M[N - 2] * tempMatrixs.combinedX.transpose())
                              .transpose()
                              .rightCols(N - 1)
                              .array()
                              .square() +
                      (params.M[N - 2] * tempMatrixs.combinedY.transpose())
                              .transpose()
                              .rightCols(N - 1)
                              .array()
                              .square() +
                      (params.M[N - 2] * tempMatrixs.combinedZ.transpose())
                              .transpose()
                              .rightCols(N - 1)
                              .array()
                              .square());

    tempMatrixs.dp =
            ((params.M_inv_t[N - 2] * tempMatrixs.combinedPX.transpose())
                     .transpose()
                     .rightCols(N - 1)
                     .array()
                     .square() +
             (params.M_inv_t[N - 2] * tempMatrixs.combinedPY.transpose())
                     .transpose()
                     .rightCols(N - 1)
                     .array()
                     .square() +
             (params.M_inv_t[N - 2] * tempMatrixs.combinedPZ.transpose())
                     .transpose()
                     .rightCols(N - 1)
                     .array()
                     .square());


    for (auto i = 0; i < N - 1; ++i) {
        tempMatrixs.diff +=
                (tempMatrixs.dr.col(i).array() / params.SigArray[i] /
                         params.SigArray[i] +
                 tempMatrixs.dp.col(i).array() * params.SigArray[i] *
                         params.SigArray[i] / hbar2)
                        .matrix();
    }
    targetParticles.col(10) =
            params.gc * params.probFactor *
            tempMatrixs.combinedProbability.rowwise().prod().array() *
            (-tempMatrixs.diff.array()).exp();
}