//
// Created by mafu on 4/4/2024.
//
#include "../headers/ClusterCalculator.h"
#include "../headers/MatrixMethod.h"
#include "../headers/ThreadPool.h"
#include "../headers/chrono.h"
#include "../headers/ctpl.h"
#include "../headers/logger.h"
#include <numeric>


std::mutex Coal::ClusterCalculator::mtx;

Eigen::MatrixXd
Coal::ClusterCalculator::Coalescence(const EventsMap &allEvents,
                                     const ClusterParams &params,
                                     std::optional<unsigned int> max_threads) {
    if (!max_threads.has_value()) {
        max_threads = 1;
    } else if (*max_threads > 0) {
        max_threads       = *max_threads - 1;
        const auto logger = Logger::getInstance().get();
        logger->info("max threads: {}", *max_threads);
    }
    // return processEventstest(allEvents, params, resolution, *max_threads);
    return processEventsV2(allEvents, params, *max_threads);
}

Coal::ClusterCalculator::MatrixMemPool::MatrixMemPool(double *memPool,
                                                      const int size_last, const int N)
    : combinedX(memPool, size_last, N), combinedY(memPool += size_last * N, size_last, N),
      combinedZ(memPool += size_last * N, size_last, N),
      combinedT(memPool += size_last * N, size_last, N),
      combinedPX(memPool += size_last * N, size_last, N),
      combinedPY(memPool += size_last * N, size_last, N),
      combinedPZ(memPool += size_last * N, size_last, N),
      combinedP0(memPool += size_last * N, size_last, N),
      combinedMass(memPool += size_last * N, size_last, N),
      combinedProbability(memPool += size_last * N, size_last, N),
      temReplicatedMaxCoeff(memPool += size_last * N, size_last, N),
      targetParticles(memPool += size_last * N, size_last, 11),
      dr(memPool += size_last * 11, size_last, N - 1),
      dp(memPool += size_last * (N - 1), size_last, N - 1),
      beta_x(memPool += size_last * (N - 1), size_last),
      beta_y(memPool += size_last, size_last), beta_z(memPool += size_last, size_last),
      diff(memPool += size_last, size_last) {}

void Coal::ClusterCalculator::MatrixMemPool::reset() {
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
    targetParticles.setZero();
    dr.setZero();
    dp.setZero();
    beta_x.setZero();
    beta_y.setZero();
    beta_z.setZero();
    diff.setZero();
}

Eigen::MatrixXd Coal::ClusterCalculator::processEventstest(
        const EventsMap &allEvents, const ClusterParams &params, ResParamsMap &resolution,
        const unsigned int num_threads) {

    constexpr long initRows = 100000;
    std::vector<long> usedRows(num_threads, 0);
    std::vector threadOutputs(num_threads, Eigen::MatrixXd(initRows, 11));
    std::vector<int> eventIDList;
    std::ranges::transform(allEvents, std::back_inserter(eventIDList),
                           [](const auto &pair) { return pair.first; });

    ThreadPool pool(num_threads);

    const SystemTimePoint start_time = SystemClock::now();

    std::atomic taskCount(0);
    for (auto i = 0; i < params.Loop; ++i) {
        auto Cell    = selectEvents(allEvents, params, resolution, i, eventIDList);
        auto subCell = selectParticles(Cell, params);
        if (subCell.empty()) {
            continue;
        }

        int taskIndex = taskCount++ % num_threads;
        pool.enqueueTask([&, subCell, taskIndex, i]() {
            processSubCellTest(subCell, params, threadOutputs[taskIndex],
                               usedRows[taskIndex], i);
        });
    }

    pool.stop();
    const long totalSize = std::accumulate(usedRows.begin(), usedRows.end(), 0L);
    Eigen::Matrix<double, Eigen::Dynamic, 11> targetParticles(totalSize, 11);
    long current_row = 0;
    for (int i = 0; i < num_threads; ++i) {
        targetParticles.block(current_row, 0, usedRows[i], 11) =
                threadOutputs[i].topRows(usedRows[i]);
        current_row += usedRows[i];
    }

    const SystemTimePoint end_time = SystemClock::now();
    const SystemTimeSpan time_span = end_time - start_time;
    printTime(time_span);
    return targetParticles;
}
void Coal::ClusterCalculator::processSubCellTest(
        const std::vector<Eigen::MatrixXd> &subCell, const ClusterParams &params,
        Eigen::MatrixXd &targetParticles, long &usedRows, int currentIndex) {

    const auto logger = Logger::getInstance().get();
    auto thread_id    = std::this_thread::get_id();
    try {
        const auto resultParticles = mainLoop(subCell, params);
        if (const long newRequiredSize = usedRows + resultParticles.rows();
            newRequiredSize > targetParticles.rows()) {
            constexpr long resize_size = 100000;
            targetParticles.conservativeResize(newRequiredSize + resize_size, 11);
        }
        targetParticles.block(usedRows, 0, resultParticles.rows(), 11) = resultParticles;
        usedRows += resultParticles.rows();
    } catch (const std::exception &e) {
        logger->error("Thread {} encountered an error at index {}: {}", thread_id,
                      currentIndex, e.what());
    }
}
Eigen::MatrixXd Coal::ClusterCalculator::processEventsV2(const EventsMap &allEvents,
                                                         const ClusterParams &params,
                                                         const unsigned int num_threads) {
    Eigen::Matrix<double, Eigen::Dynamic, 11> targetParticles(100000, 11);
    std::vector<int> eventIDList;
    std::ranges::transform(allEvents, std::back_inserter(eventIDList),
                           [](const auto &pair) { return pair.first; });
    const auto logger = Logger::getInstance().get();
    long current_row  = 0;

    const auto start_time = SystemClock::now();

    for (auto i = 0; i < params.Loop; ++i) {
        auto Cell    = selectEvents_v2(allEvents, params, eventIDList);
        auto subCell = selectParticles(Cell, params);
        if (subCell.empty()) {
            continue;
        }
        const auto result = mainLoopV1(subCell, params, num_threads);
        appendMatrixToBottom(result, targetParticles, current_row);
        logger->info("loop: {} size: {}", i, result.rows());
    }

    const auto end_time = SystemClock::now();

    printTime(end_time - start_time);
    return targetParticles.topRows(current_row);
}


auto Coal::ClusterCalculator::getjobRange(const int total_size, const int num_threads,
                                          const int thread_id) -> std::pair<int, int> {
    const int base_load = total_size / num_threads;
    const int remainder = total_size % num_threads;
    int start_index     = thread_id * base_load + std::min(thread_id, remainder);
    int end_index       = start_index + base_load + (thread_id < remainder ? 1 : 0);
    return {start_index, end_index};
}

Eigen::MatrixXd Coal::ClusterCalculator::processEvents(const EventsMap &allEvents,
                                                       const ClusterParams &params,
                                                       ResParamsMap &resolution,
                                                       const unsigned int num_threads) {
    Eigen::MatrixXd targetParticles{};

    std::vector<int> eventIDList;
    std::ranges::transform(allEvents, std::back_inserter(eventIDList),
                           [](const auto &pair) { return pair.first; });

    ThreadPool pool(num_threads);

    const SystemTimePoint start_time = SystemClock::now();

    for (auto i = 0; i < params.Loop; ++i) {
        auto Cell    = selectEvents(allEvents, params, resolution, i, eventIDList);
        auto subCell = selectParticles(Cell, params);
        if (subCell.empty()) {
            continue;
        }

        if (num_threads > 1) {
            pool.enqueueTask([subCell, &params, &targetParticles, &resolution, i]() {
                processSubCell(subCell, params, targetParticles, resolution, i);
            });
        } else {
            processSubCell(subCell, params, targetParticles, resolution, i);
        }
    }

    if (num_threads > 1) {
        pool.stop();
    }

    const SystemTimePoint end_time = SystemClock::now();
    const SystemTimeSpan time_span = end_time - start_time;
    printTime(time_span);

    return targetParticles;
}

void Coal::ClusterCalculator::processSubCell(const std::vector<Eigen::MatrixXd> &subCell,
                                             const ClusterParams &params,
                                             Eigen::MatrixXd &targetParticles,
                                             ResParamsMap &resolution,
                                             const int currentIndex) {

    const auto logger = Logger::getInstance().get();
    auto thread_id    = std::this_thread::get_id();
    try {
        const auto resultParticles = mainLoop(subCell, params);
        std::lock_guard lock(mtx);
        targetParticles.conservativeResize(
                targetParticles.rows() + resultParticles.rows(), 11);
        targetParticles.bottomRows(resultParticles.rows()) = resultParticles;

        if (params.mixEvents == 1) {
            resolution.selectEventID[currentIndex].second =
                    static_cast<int>(resultParticles.rows());
        }
    } catch (const std::exception &e) {
        logger->error("Thread {} encountered an error at index {}: {}", thread_id,
                      currentIndex, e.what());
    }
}

Eigen::MatrixXd
Coal::ClusterCalculator::mainLoop(const std::vector<Eigen::MatrixXd> &MArray,
                                  const ClusterParams &params) {
    const auto N        = params.NBody;
    const double factor = params.eventFactor;
    const int size_last = static_cast<int>(MArray[N - 1].rows());

    std::vector<bool> isValidCombination;

    if (N < 3) {
        isValidCombination.resize(1, true);
    } else
        isValidCombination.resize(N - 2, true);

    const long MAX_SIZE = (size_last * N * 11)//combined Matrix
                          + (size_last * 11)  //targetparticles
                          + (size_last * 4)   //beta_x, beta_y, beta_z, diff
                          + (size_last * (N - 1) * 2);
    std::vector<double, Eigen::aligned_allocator<double>> memPool(MAX_SIZE, 0.0);
    const auto memPtr = memPool.data();
    MatrixMemPool tempMatrixs(memPtr, size_last, N);

    const auto logger = Logger::getInstance().get();
    std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
            lorentzMatrixs(size_last, Eigen::Matrix4d::Zero());
    std::vector counts(N, 0);
    std::vector maxIndex(N, 0);
    for (auto i = 0; i < N; ++i) {
        counts[i]   = static_cast<int>(MArray[i].rows());
        maxIndex[i] = counts[i] - 1;
    }
    constexpr long init_size = 100000;
    Eigen::MatrixXd result(init_size, 11);
    long current_row = 0;
    std::unordered_map<std::string, double> distance_cache;
    double yieldAll    = 0.0;
    double yieldSelect = 0.0;
    long long loop     = 0;
    std::vector multiIndex(N, 0);
    Eigen::MatrixXd particles(N, 11);
    Eigen::MatrixXd matrixSwap(size_last, 11);

    while (true) {

        auto jumpMultiIndex = multiIndex;
        for (auto j = 0; j < N; ++j) {
            particles.row(j) = MArray[j].row(multiIndex[j]);
        }

        std::fill(isValidCombination.begin(), isValidCombination.end(), true);
        if (N > 2) {
            CheckingPortfolioValidity(particles, params, counts, multiIndex,
                                      jumpMultiIndex, distance_cache, isValidCombination);
        }

        if (std::ranges::all_of(isValidCombination, [](const bool x) { return x; })) {

            tempMatrixs.reset();
            matrixSwap.setZero();
            setMatrix(tempMatrixs, particles.topRows(N - 1), MArray[N - 1], params.pdg);
            vectorizationWithLastArray(params, tempMatrixs, lorentzMatrixs);
            yieldAll += tempMatrixs.targetParticles.col(10).sum();

            conditionSelect(tempMatrixs.targetParticles, params.probabilityCut,
                            matrixSwap);
            appendMatrixToBottom(matrixSwap, result, current_row);

            multiIndex = jumpValidLoop(multiIndex, counts, N - 3);
            if (multiIndex == maxIndex) {
                break;
            }
            loop++;
        } else {
            if (jumpMultiIndex == maxIndex) {
                break;
            }
            multiIndex = jumpMultiIndex;
        }
    }

    yieldSelect = result.block(0, 10, current_row, 1).sum();
    if (yieldSelect != 0) {
        result.block(0, 10, current_row, 1) *= yieldAll / yieldSelect;
    }
    logger->info("loop: {} yieldAll: {} yieldSelect: {} factor: {} size: {}", loop,
                 yieldAll, yieldSelect, factor, current_row);
    return result.topRows(current_row);
}
Eigen::MatrixXd
Coal::ClusterCalculator::mainLoopV1(const std::vector<Eigen::MatrixXd> &MArray,
                                    const ClusterParams &params,
                                    const unsigned int num_threads) {
    const auto size_first = static_cast<int>(MArray[0].rows());
    const int base_load   = size_first / static_cast<int>(num_threads);
    const int remainder   = size_first % static_cast<int>(num_threads);
    ThreadPool pool(num_threads);
    // ctpl::thread_pool pool(num_threads);
    constexpr long initRows = 100000;
    std::vector<long> usedRows(num_threads, 0);
    std::vector threadOutputs(num_threads, Eigen::MatrixXd(initRows, 11));
    // std::future<void> futures[num_threads];

    for (auto i = 0; i < num_threads; ++i) {
        int start_index = i * base_load + std::min(i, remainder);
        int end_index   = start_index + base_load + (i < remainder ? 1 : 0);
        pool.enqueueTask([&, start_index, end_index, i]() {
            mainLoopV2(MArray, params, start_index, end_index, threadOutputs[i],
                       usedRows[i]);
            //  pool.push([&, start_index, end_index, i](int) {
            //     mainLoopV2(MArray, params, start_index, end_index, threadOutputs[i],
            //                usedRows[i]);
        });
    }
    pool.stop();
    const long totalSize = std::accumulate(usedRows.begin(), usedRows.end(), 0L);
    Eigen::Matrix<double, Eigen::Dynamic, 11> targetParticles(totalSize, 11);
    long current_row = 0;
    for (int i = 0; i < num_threads; ++i) {
        targetParticles.block(current_row, 0, usedRows[i], 11) =
                threadOutputs[i].topRows(usedRows[i]);
        current_row += usedRows[i];
    }
    return targetParticles;
}

void Coal::ClusterCalculator::mainLoopV2(const std::vector<Eigen::MatrixXd> &MArray,
                                         const ClusterParams &params, const int start_idx,
                                         const int endidx, Eigen::MatrixXd &threadOutputs,
                                         long &usedRows) {
    const auto N        = params.NBody;
    const int size_last = static_cast<int>(MArray[N - 1].rows());
    std::vector counts(N, 0);
    std::vector maxIndex(N, 0);
    for (auto i = 0; i < N; ++i) {
        counts[i]   = static_cast<int>(MArray[i].rows());
        maxIndex[i] = counts[i] - 1;
    }
    std::vector<bool> isValidCombination;

    if (N < 3) {
        isValidCombination.resize(1, true);
    } else
        isValidCombination.resize(N - 2, true);
    const long MAX_SIZE = (size_last * N * 11)//combined Matrix
                          + (size_last * 11)  //targetparticles
                          + (size_last * 4)   //beta_x, beta_y, beta_z, diff
                          + (size_last * (N - 1) * 2);
    std::vector<double, Eigen::aligned_allocator<double>> memPool(MAX_SIZE, 0.0);
    const auto memPtr = memPool.data();
    MatrixMemPool tempMatrixs(memPtr, size_last, N);

    std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
            lorentzMatrixs(size_last, Eigen::Matrix4d::Zero());
    std::unordered_map<std::string, double> distance_cache;
    double yieldAll = 0.0;
    long long loop  = 0;
    std::vector multiIndex(N, 0);
    multiIndex[0] = start_idx;
    maxIndex[0]   = endidx - 1;
    counts[0]     = endidx;
    Eigen::MatrixXd particles(N, 11);
    Eigen::MatrixXd matrixSwap(size_last, 11);
    while (true) {
        auto jumpMultiIndex = multiIndex;

        for (auto j = 0; j < N; ++j) {
            particles.row(j) = MArray[j].row(multiIndex[j]);
        }

        std::fill(isValidCombination.begin(), isValidCombination.end(), true);
        if (N > 2) {
            CheckingPortfolioValidity(particles, params, counts, multiIndex,
                                      jumpMultiIndex, distance_cache, isValidCombination);
        }

        if (std::ranges::all_of(isValidCombination, [](const bool x) { return x; })) {
            tempMatrixs.reset();
            matrixSwap.setZero();
            setMatrix(tempMatrixs, particles.topRows(N - 1), MArray[N - 1], params.pdg);
            vectorizationWithLastArray(params, tempMatrixs, lorentzMatrixs);
            yieldAll += tempMatrixs.targetParticles.col(10).sum();
            conditionSelect(tempMatrixs.targetParticles, params.probabilityCut,
                            matrixSwap);
            // select(tempMatrixs.targetParticles, 10, params.probabilityCut, matrixSwap);
            appendMatrixToBottom(matrixSwap, threadOutputs, usedRows);
            multiIndex = jumpValidLoop(multiIndex, counts, N - 3);
            if (multiIndex == maxIndex) {
                break;
            }
            loop++;
        } else {
            if (jumpMultiIndex == maxIndex) {
                break;
            }
            multiIndex = jumpMultiIndex;
        }
    }
    if (const double yieldSelect = threadOutputs.block(0, 10, usedRows, 1).sum();
        yieldSelect != 0) {
        threadOutputs.block(0, 10, usedRows, 1) *= yieldAll / yieldSelect;
    }
}

std::vector<int>
Coal::ClusterCalculator::jumpValidLoop(const std::vector<int> &multiIndex,
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

    if (std::ranges::all_of(newMultiIndex, [](const int x) { return x == 0; })) {
        for (int i = 0; i < newMultiIndex.size(); ++i) {
            newMultiIndex[i] = counts[i] - 1;
        }
    }

    return newMultiIndex;
}

std::string
Coal::ClusterCalculator::createKeyFromMultiIndex(const std::vector<int> &multiIndex,
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

void Coal::ClusterCalculator::CheckingPortfolioValidity(
        const Eigen::MatrixXd &ParticlesList, const ClusterParams &params,
        const std::vector<int> &counts, const std::vector<int> &multIndex,
        std::vector<int> &jumpMultiIndex,
        std::unordered_map<std::string, double> &distanceCache,
        std::vector<bool> &isValidCombination) {

    const auto N = static_cast<int>(ParticlesList.rows());

    for (auto i = 1; i < N - 1; ++i) {
        auto dis_temp   = 0.0;
        std::string key = createKeyFromMultiIndex(multIndex, i);

        if (auto it = distanceCache.find(key); it != distanceCache.end()) {
            dis_temp = it->second;
        } else {
            for (auto index = distanceCache.begin(); index != distanceCache.end();) {
                if (std::ranges::count(index->first, '-') == i) {
                    index = distanceCache.erase(index);
                } else {
                    ++index;
                }
            }

            Eigen::Matrix<double, Eigen::Dynamic, 11> tempParticle =
                    ParticlesList.topRows(i + 1);
            boostToComMatrix(tempParticle);

            auto [diff_r, diff_p] = JacobiCoordinatesMatrix_Vec(tempParticle, params);

            for (auto j = 0; j < i; ++j) {
                dis_temp += (diff_r(j) * diff_r(j) * params.Sig_inv_sq[j] +
                             diff_p(j) * diff_p(j) * params.Sig_sq_hbar2[j]);
            }
            distanceCache[key] = dis_temp;
        }

        if (dis_temp > (params.precision + (i - 1) * 8)) {
            isValidCombination[i - 1] = false;
            jumpMultiIndex            = jumpValidLoop(multIndex, counts, i - 1);
            return;
        }
    }
}
void Coal::ClusterCalculator::setMatrix(MatrixMemPool &temMatrix,
                                        const Eigen::MatrixXd &particles,
                                        const Eigen::MatrixXd &lastParticles,
                                        const int pdg) {
    const int size      = static_cast<int>(particles.rows());
    const int size_last = static_cast<int>(lastParticles.rows());
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

    temMatrix.targetParticles.col(0) = pdg * Eigen::VectorXd::Ones(size_last);
    temMatrix.targetParticles.col(1) = temMatrix.combinedPX.rowwise().sum();
    temMatrix.targetParticles.col(2) = temMatrix.combinedPY.rowwise().sum();
    temMatrix.targetParticles.col(3) = temMatrix.combinedPZ.rowwise().sum();
    temMatrix.targetParticles.col(4) = temMatrix.combinedP0.rowwise().sum();
    temMatrix.targetParticles.col(5) = temMatrix.combinedMass.rowwise().sum();
    temMatrix.targetParticles.col(6) = temMatrix.combinedX.rowwise().mean();
    temMatrix.targetParticles.col(7) = temMatrix.combinedY.rowwise().mean();
    temMatrix.targetParticles.col(8) = temMatrix.combinedZ.rowwise().mean();
    temMatrix.targetParticles.col(9) = temMatrix.combinedT.rowwise().maxCoeff();
}

void Coal::ClusterCalculator::vectorizationWithLastArray(
        const ClusterParams &params, MatrixMemPool &tempMatrixs,
        std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
                &lorentzMatrixs) {

    const auto N = params.NBody;
    for (auto &matrix: lorentzMatrixs) {
        matrix.setZero();
    }

    tempMatrixs.beta_x = tempMatrixs.targetParticles.col(1).array() /
                         tempMatrixs.targetParticles.col(4).array();
    tempMatrixs.beta_y = tempMatrixs.targetParticles.col(2).array() /
                         tempMatrixs.targetParticles.col(4).array();
    tempMatrixs.beta_z = tempMatrixs.targetParticles.col(3).array() /
                         tempMatrixs.targetParticles.col(4).array();

    calculateLorentz(tempMatrixs.beta_x, tempMatrixs.beta_y, tempMatrixs.beta_z,
                     lorentzMatrixs);

    applyLorentzBoost(tempMatrixs.combinedX, tempMatrixs.combinedY, tempMatrixs.combinedZ,
                      tempMatrixs.combinedT, tempMatrixs.combinedPX,
                      tempMatrixs.combinedPY, tempMatrixs.combinedPZ,
                      tempMatrixs.combinedP0, lorentzMatrixs);
    tempMatrixs.temReplicatedMaxCoeff =
            tempMatrixs.combinedT.rowwise()
                    .maxCoeff()
                    .replicate(1, tempMatrixs.combinedT.cols())
                    .array() -
            tempMatrixs.combinedT.array();

    const auto commonFactor =
            tempMatrixs.temReplicatedMaxCoeff.array() / tempMatrixs.combinedP0.array();

    // tempMatrixs.combinedX = tempMatrixs.temReplicatedMaxCoeff.array() *
    //                                 tempMatrixs.combinedPX.array() /
    //                                 tempMatrixs.combinedP0.array() +
    //                         tempMatrixs.combinedX.array();
    // tempMatrixs.combinedY = tempMatrixs.temReplicatedMaxCoeff.array() *
    //                                 tempMatrixs.combinedPY.array() /
    //                                 tempMatrixs.combinedP0.array() +
    //                         tempMatrixs.combinedY.array();
    // tempMatrixs.combinedZ = tempMatrixs.temReplicatedMaxCoeff.array() *
    //                                 tempMatrixs.combinedPZ.array() /
    //                                 tempMatrixs.combinedP0.array() +
    //                         tempMatrixs.combinedZ.array();
    tempMatrixs.combinedX.array() += commonFactor * tempMatrixs.combinedPX.array();
    tempMatrixs.combinedY.array() += commonFactor * tempMatrixs.combinedPY.array();
    tempMatrixs.combinedZ.array() += commonFactor * tempMatrixs.combinedPZ.array();

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

    tempMatrixs.dp = ((params.M_inv_t[N - 2] * tempMatrixs.combinedPX.transpose())
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
        tempMatrixs.diff += (tempMatrixs.dr.col(i).array() * params.Sig_inv_sq[i] +
                             tempMatrixs.dp.col(i).array() * params.Sig_sq_hbar2[i])
                                    .matrix();
    }
    tempMatrixs.targetParticles.col(10) =
            params.gc * params.probFactor *
            tempMatrixs.combinedProbability.rowwise().prod().array() *
            (-tempMatrixs.diff.array()).exp();
}
