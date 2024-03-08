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

Coal::ComputationMatrices::ComputationMatrices(double *memPool,
                                               const int size_last, const int N)
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
    std::map<std::pair<double, double>, std::vector<double>> ptArray;
    std::map<std::pair<double, double>, double> yeildArray;
    const RapidityArray rapidityArray =
            output["RapidityRange"].as<RapidityArray>();
    const bool extended = output["Extended"].as<bool>();
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(params.ptBins.second, 0.0);
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
                               yeildArray, params, rapidityArray, extended);

    for (const auto &rap: rapidityArray) {
        ptOutput << "Rapidity range: " << rap.first << "<y<" << rap.second
                 << ", cluster yield:" << yeildArray[rap] << "\n";
        for (auto i = 0; i < params.ptBins.second; ++i) {
            double pt = params.ptBins.first / 2 +
                        static_cast<double>(i) * params.ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * params.ptBins.first *
                                std::abs((rap.second - rap.first)));
            ptOutput << pt << " " << ptArray[rap][i] << "\n";
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
void Coal::clusterNoThread(const EventsMap &allEvents,
                           const std::string &outputFilename,
                           const std::string &ptFilename,
                           const ClusterParams &params,
                           const RapidityArray &rapidityArray) {
    std::ofstream clusterOutput(outputFilename);
    std::ofstream ptOutput(ptFilename);
    std::map<std::pair<double, double>, std::vector<double>> ptArray;
    std::map<std::pair<double, double>, double> yeildArray;
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(params.ptBins.second, 0.0);
        yeildArray[rap] = 0.0;
    }

    ParticleArray targetParticles{};
    ParticleArray resultParticles{};
    auto start_time = std::chrono::steady_clock::now();

    for (auto i = 0; i < params.Loop; ++i) {
        auto Cell    = selectEvents(allEvents, params);
        auto subCell = selectParticles(Cell, params);
        if (subCell.empty()) {
            continue;
        }

        resultParticles = mainLoopV2(subCell, params);
        targetParticles.insert(targetParticles.end(),
                               std::make_move_iterator(resultParticles.begin()),
                               std::make_move_iterator(resultParticles.end()));
        resultParticles.clear();
    }

    outputTargetParticle(targetParticles, clusterOutput, ptArray, yeildArray,
                         params, rapidityArray);
    for (const auto &rap: rapidityArray) {
        ptOutput << "Rapidity range: " << rap.first << "<y<" << rap.second
                 << ", cluster yield:" << yeildArray[rap] << "\n";
        for (auto i = 0; i < params.ptBins.second; ++i) {
            double pt = params.ptBins.first / 2 +
                        static_cast<double>(i) * params.ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * params.ptBins.first *
                                std::abs((rap.second - rap.first)));
            ptOutput << pt << " " << ptArray[rap][i] << "\n";
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
void Coal::clusterMix(const EventsMap &allEvents,
                      const std::string &outputFilename,
                      const std::string &ptFilename,
                      const ClusterParams &params,
                      const RapidityArray &rapidityArray) {
    std::ofstream clusterOutput(outputFilename);
    std::ofstream ptOutput(ptFilename);
    std::map<std::pair<double, double>, std::vector<double>> ptArray;
    std::map<std::pair<double, double>, double> yeildArray;
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(params.ptBins.second, 0.0);
        yeildArray[rap] = 0.0;
    }

    ParticleArray targetParticles{};
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (max_threads > 0)
        max_threads -= 1;

    std::vector<std::thread> threads;
    std::queue<MultiParticleArray> taskQueue{};

    auto start_time = std::chrono::steady_clock::now();

    for (auto i = 0; i < params.Loop; ++i) {
        auto Cell    = selectEvents(allEvents, params);
        auto subCell = selectParticles(Cell, params);
        if (subCell.empty()) {
            continue;
        }
        taskQueue.push(subCell);
    }
    while (!taskQueue.empty()) {
        while (threads.size() < max_threads && !taskQueue.empty()) {
            auto subCell = taskQueue.front();
            taskQueue.pop();
            threads.emplace_back(processSubCell, subCell, params,
                                 std::ref(targetParticles));
        }
        for (auto &thread: threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        std::erase_if(threads,
                      [](const std::thread &x) { return !x.joinable(); });
    }

    outputTargetParticle(targetParticles, clusterOutput, ptArray, yeildArray,
                         params, rapidityArray);
    for (const auto &rap: rapidityArray) {
        ptOutput << "Rapidity range: " << rap.first << "<y<" << rap.second
                 << ", cluster yield:" << yeildArray[rap] << "\n";
        for (auto i = 0; i < params.ptBins.second; ++i) {
            double pt = params.ptBins.first / 2 +
                        static_cast<double>(i) * params.ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * params.ptBins.first *
                                std::abs((rap.second - rap.first)));
            ptOutput << pt << " " << ptArray[rap][i] << "\n";
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
void Coal::clusterThreadMatrix(const EventsMap &allEvents,
                               const std::string &outputFilename,
                               const std::string &ptFilename,
                               const ClusterParams &params,
                               const YAML::Node &output) {
    std::ofstream clusterOutput(outputFilename);
    std::ofstream ptOutput(ptFilename);
    std::map<std::pair<double, double>, std::vector<double>> ptArray;
    std::map<std::pair<double, double>, double> yeildArray;
    const RapidityArray rapidityArray =
            output["RapidityRange"].as<RapidityArray>();
    const bool extended = output["Extended"].as<bool>();
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(params.ptBins.second, 0.0);
        yeildArray[rap] = 0.0;
    }

    Eigen::MatrixXd targetParticles{};
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (max_threads > 0)
        max_threads -= 1;

    std::vector<std::thread> threads;
    std::queue<std::vector<Eigen::MatrixXd>> taskQueue{};

    auto start_time = std::chrono::steady_clock::now();

    for (auto i = 0; i < params.Loop; ++i) {
        auto Cell    = selectEvents(allEvents, params);
        auto subCell = selectParticlesMatrix(Cell, params);
        if (subCell.empty()) {
            continue;
        }
        taskQueue.push(subCell);
    }
    while (!taskQueue.empty()) {
        while (threads.size() < max_threads && !taskQueue.empty()) {
            auto subCell = taskQueue.front();
            taskQueue.pop();
            threads.emplace_back(processSubCellMatrix, subCell, params,
                                 std::ref(targetParticles));
        }
        for (auto &thread: threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        std::erase_if(threads,
                      [](const std::thread &x) { return !x.joinable(); });
    }

    outputTargetParticleMatrix(targetParticles, clusterOutput, ptArray,
                               yeildArray, params, rapidityArray, extended);
    for (const auto &rap: rapidityArray) {
        ptOutput << "Rapidity range: " << rap.first << "<y<" << rap.second
                 << ", cluster yield:" << yeildArray[rap] << "\n";
        for (auto i = 0; i < params.ptBins.second; ++i) {
            double pt = params.ptBins.first / 2 +
                        static_cast<double>(i) * params.ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * params.ptBins.first *
                                std::abs((rap.second - rap.first)));
            ptOutput << pt << " " << ptArray[rap][i] << "\n";
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
    std::unique_lock lock(queue_mutex);
    while (!taskQueue.empty()) {
        taskQueue.pop();
    }
    stop_all = false;
    lock.unlock();

    cv.notify_all();// 确保所有等待的线程都被唤醒

    for (auto &worker: workers) {
        if (worker.joinable()) {
            worker.join();// 等待所有旧线程结束
        }
    }
    workers.clear();// 清空线程向量以便重新创建线程
}
void Coal::clusterThreadV2(const EventsMap &allEvents,
                           const std::string &outputFilename,
                           const std::string &ptFilename,
                           const ClusterParams &params,
                           const YAML::Node &config) {
    resetGlobalState();
    std::ofstream clusterOutput(outputFilename);
    std::ofstream ptOutput(ptFilename);
    std::map<std::pair<double, double>, std::vector<double>> ptArray;
    std::map<std::pair<double, double>, double> yeildArray;
    auto output = config["Output"];
    const RapidityArray rapidityArray =
            output["RapidityRange"].as<RapidityArray>();
    const bool extended = output["Extended"].as<bool>();
    for (const auto &rap: rapidityArray) {
        ptArray[rap]    = std::vector(params.ptBins.second, 0.0);
        yeildArray[rap] = 0.0;
    }

    Eigen::MatrixXd targetParticles{};
    unsigned int max_threads =
            config["General"]["Parallel"]["Cores"].as<unsigned>();
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
                               yeildArray, params, rapidityArray, extended);
    for (const auto &rap: rapidityArray) {
        ptOutput << "Rapidity range: " << rap.first << "<y<" << rap.second
                 << ", cluster yield:" << yeildArray[rap] << "\n";
        for (auto i = 0; i < params.ptBins.second; ++i) {
            double pt = params.ptBins.first / 2 +
                        static_cast<double>(i) * params.ptBins.first;
            ptArray[rap][i] /= (2 * M_PI * pt * params.ptBins.first *
                                std::abs((rap.second - rap.first)));
            ptOutput << pt << " " << ptArray[rap][i] << "\n";
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

void Coal::processSubCell(const MultiParticleArray &subCell,
                          const ClusterParams &params,
                          ParticleArray &targetParticles) {
    auto resultParticles = mainLoopV2(subCell, params);
    std::lock_guard lock(mtx);
    targetParticles.insert(targetParticles.end(), resultParticles.begin(),
                           resultParticles.end());
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
        std::map<std::pair<double, double>, std::vector<double>> &ptArray,
        std::map<std::pair<double, double>, double> &yeildArray,
        const ClusterParams &params, const RapidityArray &rapidityArray,
        const bool &extended) {
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
            if (rapidity > rap.first && rapidity < rap.second) {
                if (const int index = static_cast<int>(pt / d_pt);
                    index < ptBins) {
                    ptArray[rap][index] += targetParticles(i, 10) / totalEvents;
                }
                yeildArray[rap] += targetParticles(i, 10) / totalEvents;
            }
        }
    }
    if (extended) {
        outputCluster << targetParticles << "\n";
    }
    outputCluster.close();
}
void Coal::outputTargetParticle(
        const ParticleArray &targetParticles, std::ofstream &outputCluster,
        std::map<std::pair<double, double>, std::vector<double>> &ptArray,
        std::map<std::pair<double, double>, double> &yeildArray,
        const ClusterParams &params, const RapidityArray &rapidityArray) {
    const double d_pt        = params.ptBins.first;
    const int ptBins         = params.ptBins.second;
    const double totalEvents = params.eventFactor;
    outputCluster << "#events:" << totalEvents << " "
                  << "size:" << targetParticles.size() << "\n";
    for (auto &cluster: targetParticles) {
        const double pt       = cluster.getPT();
        const double rapidity = cluster.getRapidity();
        for (const auto &rap: rapidityArray) {
            if (rapidity > rap.first && rapidity < rap.second) {
                if (const int index = static_cast<int>(pt / d_pt);
                    index < ptBins) {
                    ptArray[rap][index] += cluster.probability / totalEvents;
                }
                yeildArray[rap] += cluster.probability / totalEvents;
            }
        }
        outputCluster << params.pdg << " " << cluster.px << " " << cluster.py
                      << " " << cluster.pz << " " << cluster.mass << " "
                      << cluster.x << " " << cluster.y << " " << cluster.z
                      << " " << cluster.freeze_out_time << " " << cluster.p0
                      << " " << cluster.probability << "\n";
    }
    outputCluster.close();
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
Coal::ParticleArray
Coal::mainLoopV2(const MultiParticleArray &ReactionParticles,
                 const ClusterParams &params) {
    const auto N        = params.NBody;
    const double factor = params.eventFactor;
    std::vector counts(N, 0);
    std::vector maxIndex(N, 0);
    for (auto i = 0; i < N; ++i) {
        counts[i]   = static_cast<int>(ReactionParticles[i].size());
        maxIndex[i] = counts[i] - 1;
    }

    ParticleArray targetParticleArray{};
    std::unordered_map<std::string, double> distance_cache;
    double yieldAll    = 0.0;
    double yieldSelect = 0.0;
    long long loop     = 0;
    std::vector multiIndex(N, 0);
    while (true) {
        auto jumpMultiIndex = multiIndex;
        ParticleArray particles{};
        for (auto j = 0; j < N; ++j) {
            particles.push_back(ReactionParticles[j][multiIndex[j]]);
        }
        Particle targetParticle{};

        nBodyCoalV2(particles, targetParticle, params, counts, multiIndex,
                    jumpMultiIndex, distance_cache);

        if (jumpMultiIndex != multiIndex) {
            if (jumpMultiIndex == maxIndex) {
                break;
            }
            multiIndex = jumpMultiIndex;
            continue;
        }
        yieldAll += targetParticle.probability;
        if (targetParticle.probability > 0.0 &&
            targetParticle.probability > params.probabilityCut) {
            yieldSelect += targetParticle.probability;
            targetParticleArray.push_back(targetParticle);
        }
        loop++;
        if (!incrementIndex(multiIndex, counts) || multiIndex == maxIndex) {
            break;
        }
    }
    for (auto &particle: targetParticleArray) {
        particle.probability *= yieldAll / yieldSelect;
    }
    std::cout << "loop:" << loop << " "
              << "yieldAll:" << yieldAll << " "
              << "yieldSelect:" << yieldSelect << " "
              << "factor:" << factor << " "
              << "size:" << targetParticleArray.size() << "\n";

    return targetParticleArray;
}

// Eigen::MatrixXd Coal::mainLoopMatrix(const std::vector<Eigen::MatrixXd> &MArray,
//                                      const ClusterParams &params) {
//     const auto N         = params.NBody;
//     const double factor  = params.eventFactor;
//     const auto size_last = MArray[N - 1].rows();
//     std::vector counts(N, 0);
//     std::vector maxIndex(N, 0);
//     for (auto i = 0; i < N; ++i) {
//         counts[i]   = static_cast<int>(MArray[i].rows());
//         maxIndex[i] = counts[i] - 1;
//     }
//     Eigen::MatrixXd targetParticleArray;
//     std::unordered_map<std::string, double> distance_cache;
//     double yieldAll    = 0.0;
//     double yieldSelect = 0.0;
//     long long loop     = 0;
//     std::vector multiIndex(N, 0);
//     Eigen::Matrix<double, Eigen::Dynamic, 11> targetParticle(size_last, 11);
//     while (true) {
//         auto jumpMultiIndex = multiIndex;
//         Eigen::MatrixXd particles(N, 11);
//         for (auto j = 0; j < N; ++j) {
//             particles.row(j) = MArray[j].row(multiIndex[j]);
//         }
//
//         auto isValidCombination =
//                 checkCombinationList(particles, params, counts, multiIndex,
//                                      jumpMultiIndex, distance_cache);
//         if (std::ranges::all_of(isValidCombination,
//                                 [](const bool x) { return x; })) {
//
//             targetParticle.setZero();
//             batchProcessLastParticlesCols(
//                     particles.block(0, 0, particles.rows() - 1,
//                                     particles.cols()),
//                     MArray[N - 1], params, targetParticle);
//             yieldAll += targetParticle.col(10).sum();
//             if (Eigen::Array<bool, Eigen::Dynamic, 1> mask =
//                         targetParticle.col(10).array() > params.probabilityCut;
//                 mask.any()) {
//                 std::vector<int> indices;
//                 for (int i = 0; i < mask.size(); ++i) {
//                     if (mask(i)) {
//                         indices.push_back(i);
//                     }
//                 }
//
//                 Eigen::MatrixXd filteredTargetParticle(indices.size(),
//                                                        targetParticle.cols());
//                 for (size_t i = 0; i < indices.size(); ++i) {
//                     filteredTargetParticle.row(i) =
//                             targetParticle.row(indices[i]);
//                 }
//
//                 yieldSelect += filteredTargetParticle.col(10).sum();
//
//                 targetParticleArray.conservativeResize(
//                         targetParticleArray.rows() +
//                                 filteredTargetParticle.rows(),
//                         11);
//                 targetParticleArray.bottomRows(filteredTargetParticle.rows()) =
//                         filteredTargetParticle;
//             }
//
//             multiIndex = jumpValidLoop(multiIndex, counts, N - 3);
//             if (multiIndex == maxIndex) {
//                 break;
//             }
//
//         } else {
//             if (jumpMultiIndex == maxIndex) {
//                 break;
//             }
//             multiIndex = jumpMultiIndex;
//             continue;
//         }
//         loop++;
//     }
//
//     targetParticleArray.col(10) *= yieldAll / yieldSelect;
//     std::cout << "loop:" << loop << " "
//               << "yieldAll:" << yieldAll << " "
//               << "yieldSelect:" << yieldSelect << " "
//               << "factor:" << factor << " "
//               << "size:" << targetParticleArray.rows() << "\n";
//     return targetParticleArray;
// }

Eigen::MatrixXd Coal::mainLoopMatrix(const std::vector<Eigen::MatrixXd> &MArray,
                                     const ClusterParams &params) {
    const auto N         = params.NBody;
    const double factor  = params.eventFactor;
    const auto size_last = MArray[N - 1].rows();
    const int MAX_SIZE   = (size_last * N * 10)//combined Matrix
                         + (size_last * 4)     //beta_x, beta_y, beta_z, diff
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
                for (size_t i = 0; i < indices.size(); ++i) {
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


Coal::ParticleArray Coal::mainLoop(const MultiParticleArray &ReactionParticles,
                                   const ClusterParams &params) {
    const auto N        = params.NBody;
    const double factor = params.eventFactor;
    std::vector counts(N, 0);
    long long loopMax = 1;
    for (auto i = 0; i < N; ++i) {
        counts[i] = static_cast<int>(ReactionParticles[i].size());
        loopMax *= counts[i];
    }
    ParticleArray targetParticleArray{};
    double yieldAll    = 0.0;
    double yieldSelect = 0.0;
    long long loop     = 0;
    std::vector lastmultiIndex(N, -1);
    for (long long i = 0; i < loopMax;) {

        const auto multiIndex = indexToMultiIndex(i, counts);
        auto jumpMultiIndex   = multiIndex;
        ParticleArray particles{};
        for (auto j = 0; j < N; ++j) {
            particles.push_back(ReactionParticles[j][multiIndex[j]]);
        }
        Particle targetParticle{};

        nBodyCoal(particles, targetParticle, params, counts, multiIndex,
                  jumpMultiIndex, lastmultiIndex);

        if (jumpMultiIndex != multiIndex) {
            i = multiIndexToIndex(jumpMultiIndex, counts);
            if (i >= loopMax) {
                break;
            }
            continue;
        }
        yieldAll += targetParticle.probability;

        if (targetParticle.probability > 0.0 &&
            targetParticle.probability > params.probabilityCut) {
            yieldSelect += targetParticle.probability;
            targetParticleArray.push_back(targetParticle);
        }
        loop++;
        ++i;
    }
    for (auto &particle: targetParticleArray) {
        particle.probability *= yieldAll / yieldSelect;
    }
    std::cout << "loopMax:" << loopMax << " "
              << "loopCut:" << loop << " "
              << "yieldAll:" << yieldAll << " "
              << "yieldSelect:" << yieldSelect << " "
              << "factor:" << factor << " "
              << "size:" << targetParticleArray.size() << "\n";

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


void Coal::nBodyCoal(const ParticleArray &Particles, Particle &targetParticle,
                     const ClusterParams &params,
                     const std::vector<int> &counts,
                     const std::vector<int> &multIndex,
                     std::vector<int> &jumpMultiIndex,
                     std::vector<int> &lastMultiIndex) {
    constexpr double hbar = 0.19733;

    constexpr double hbar2 = hbar * hbar;
    const auto N           = params.NBody;

    targetParticle.getResultParticleData(Particles);
    targetParticle.pdg         = params.pdg;
    const auto boost_particles = boostToCOM(Particles, targetParticle);

    auto [diff_r, diff_p] = JacobiCoordinates(boost_particles, params);

    auto dis                  = std::vector(N - 1, 0.0);
    double distotal           = 0.0;
    const bool indexesChanged = !std::equal(
            lastMultiIndex.begin(),
            lastMultiIndex.begin() + N - params.precision, multIndex.begin());

    for (auto i = 0; i < N - 1; ++i) {

        dis[i] = (diff_r[i] * diff_r[i] / params.SigArray[i] /
                          params.SigArray[i] +
                  diff_p[i] * diff_p[i] * params.SigArray[i] *
                          params.SigArray[i] / hbar2);
        if (indexesChanged && i != (N - params.precision) &&
            dis[i] > (25 + i * 8)) {
            targetParticle.probability = 0.0;
            jumpMultiIndex             = jumpValidLoop(multIndex, counts, i);
            return;
        }
        distotal += dis[i];
    }
    const double prob_factor = std::accumulate(
            boost_particles.begin(), boost_particles.end(), 1.0,
            [](const auto &a, const auto &b) { return a * b.probability; });
    targetParticle.probability =
            prob_factor * params.gc * params.probFactor * std::exp(-distotal);
}
void Coal::nBodyCoalV2(const ParticleArray &Particles, Particle &targetParticle,
                       const ClusterParams &params,
                       const std::vector<int> &counts,
                       const std::vector<int> &multIndex,
                       std::vector<int> &jumpMultiIndex,
                       std::unordered_map<std::string, double> &distanceCache) {
    constexpr double hbar  = 0.19733;
    constexpr double hbar2 = hbar * hbar;
    const auto N           = params.NBody;
    auto dis_temp          = 0.0;

    for (auto i = 1; i < N; ++i) {
        dis_temp        = 0.0;
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
            std::vector temParticles(Particles.begin(),
                                     Particles.begin() + i + 1);
            targetParticle.getResultParticleData(temParticles);
            const auto boost_particles =
                    boostToCOM(temParticles, targetParticle);

            const auto [diff_r, diff_p] =
                    JacobiCoordinates(boost_particles, params);

            for (auto j = 0; j < i; ++j) {
                dis_temp += (diff_r[j] * diff_r[j] / params.SigArray[j] /
                                     params.SigArray[j] +
                             diff_p[j] * diff_p[j] * params.SigArray[j] *
                                     params.SigArray[j] / hbar2);
            }
            distanceCache[key] = dis_temp;
        }
        if (dis_temp > (25 + (i - 1) * 8)) {
            targetParticle.probability = 0.0;
            jumpMultiIndex = jumpValidLoop(multIndex, counts, i - 1);
            return;
        }
    }
    double prob_factor;
    if (N == 2) {

        prob_factor = std::accumulate(
                Particles.begin(), Particles.end(), 1.0,
                [](const auto &a, const auto &b) { return a * b.probability; });
    } else {
        prob_factor = 1;
    }
    targetParticle.pdg = params.pdg;
    targetParticle.probability =
            prob_factor * params.gc * params.probFactor * std::exp(-dis_temp);
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
void Coal::batchProcessLastParticles(
        const Eigen::MatrixXd &Particles, const Eigen::MatrixXd &LastParticles,
        const ClusterParams &params,
        Eigen::Matrix<double, -1, 11> &targetParticles) {
    constexpr double hbar  = 0.19733;
    constexpr double hbar2 = hbar * hbar;
    const auto size        = static_cast<int>(Particles.rows());
    const auto size_last   = static_cast<int>(LastParticles.rows());
    const auto N           = params.NBody;

    Eigen::MatrixXd combinedX           = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedY           = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedZ           = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedT           = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedPX          = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedPY          = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedPZ          = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedP0          = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedMass        = Eigen::MatrixXd::Zero(size_last, N);
    Eigen::MatrixXd combinedProbability = Eigen::MatrixXd::Zero(size_last, N);

    for (auto i = 0; i < size; ++i) {
        combinedX.col(i).setConstant(Particles(i, 6));

        combinedY.col(i).setConstant(Particles(i, 7));

        combinedZ.col(i).setConstant(Particles(i, 8));
        combinedT.col(i).setConstant(Particles(i, 9));

        combinedPX.col(i).setConstant(Particles(i, 1));

        combinedPY.col(i).setConstant(Particles(i, 2));

        combinedPZ.col(i).setConstant(Particles(i, 3));

        combinedP0.col(i).setConstant(Particles(i, 5));

        combinedMass.col(i).setConstant(Particles(i, 4));

        combinedProbability.col(i).setConstant(Particles(i, 10));
    }
    combinedX.col(N - 1)           = LastParticles.col(6);
    combinedY.col(N - 1)           = LastParticles.col(7);
    combinedZ.col(N - 1)           = LastParticles.col(8);
    combinedT.col(N - 1)           = LastParticles.col(9);
    combinedPX.col(N - 1)          = LastParticles.col(1);
    combinedPY.col(N - 1)          = LastParticles.col(2);
    combinedPZ.col(N - 1)          = LastParticles.col(3);
    combinedP0.col(N - 1)          = LastParticles.col(5);
    combinedMass.col(N - 1)        = LastParticles.col(4);
    combinedProbability.col(N - 1) = LastParticles.col(10);


    Eigen::VectorXd px_total = combinedPX.rowwise().sum();
    Eigen::VectorXd py_total = combinedPY.rowwise().sum();
    Eigen::VectorXd pz_total = combinedPZ.rowwise().sum();
    Eigen::VectorXd p0_total = combinedP0.rowwise().sum();

    targetParticles.col(0) = params.pdg * Eigen::VectorXd::Ones(size_last);
    targetParticles.col(1) = px_total;
    targetParticles.col(2) = py_total;
    targetParticles.col(3) = pz_total;
    targetParticles.col(4) = p0_total;
    targetParticles.col(5) = combinedMass.rowwise().sum();
    targetParticles.col(6) = combinedX.rowwise().mean();
    targetParticles.col(7) = combinedY.rowwise().mean();
    targetParticles.col(8) = combinedZ.rowwise().mean();
    targetParticles.col(9) = combinedT.rowwise().maxCoeff();


    const Eigen::VectorXd beta_x = px_total.array() / p0_total.array();
    const Eigen::VectorXd beta_y = py_total.array() / p0_total.array();
    const Eigen::VectorXd beta_z = pz_total.array() / p0_total.array();

    const auto lambda = calculateLorentz(beta_x, beta_y, beta_z);
    applyLorentzBoost(combinedX, combinedY, combinedZ, combinedT, combinedPX,
                      combinedPY, combinedPZ, combinedP0, lambda);

    Eigen::VectorXd t_max            = combinedT.rowwise().maxCoeff();
    Eigen::MatrixXd t_max_replicated = t_max.replicate(1, combinedT.cols());
    combinedX = (t_max_replicated.array() - combinedT.array()) *
                        combinedPX.array() / combinedP0.array() +
                combinedX.array();
    combinedY = (t_max_replicated.array() - combinedT.array()) *
                        combinedPY.array() / combinedP0.array() +
                combinedY.array();
    combinedZ = (t_max_replicated.array() - combinedT.array()) *
                        combinedPZ.array() / combinedP0.array() +
                combinedZ.array();

    auto transformCoordinates = [&](const Eigen::MatrixXd &mat,
                                    const Eigen::MatrixXd &transformat) {
        return (transformat * mat.transpose())
                .transpose()
                .rightCols(mat.cols() - 1);
    };

    Eigen::MatrixXd x_jacobi = transformCoordinates(combinedX, params.M[N - 2]);
    Eigen::MatrixXd y_jacobi = transformCoordinates(combinedY, params.M[N - 2]);
    Eigen::MatrixXd z_jacobi = transformCoordinates(combinedZ, params.M[N - 2]);
    Eigen::MatrixXd px_jacobi =
            transformCoordinates(combinedPX, params.M_inv_t[N - 2]);
    Eigen::MatrixXd py_jacobi =
            transformCoordinates(combinedPY, params.M_inv_t[N - 2]);
    Eigen::MatrixXd pz_jacobi =
            transformCoordinates(combinedPZ, params.M_inv_t[N - 2]);

    Eigen::MatrixXd d_r = Eigen::MatrixXd::Zero(size_last, N - 1);
    Eigen::MatrixXd d_p = Eigen::MatrixXd::Zero(size_last, N - 1);
    d_r = x_jacobi.array().square() + y_jacobi.array().square() +
          z_jacobi.array().square();
    d_p = px_jacobi.array().square() + py_jacobi.array().square() +
          pz_jacobi.array().square();

    Eigen::VectorXd diff = Eigen::VectorXd::Zero(size_last);

    for (auto i = 0; i < N - 1; ++i) {
        diff += (d_r.col(i).array() / params.SigArray[i] / params.SigArray[i] +
                 d_p.col(i).array() * params.SigArray[i] * params.SigArray[i] /
                         hbar2)
                        .matrix();
    }
    Eigen::VectorXd probability = params.gc * params.probFactor *
                                  combinedProbability.rowwise().prod().array() *
                                  (-diff.array()).exp();

    targetParticles.col(10) = probability;
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

