//
// Created by mafu on 1/9/2024.
//
#include <numeric>

#include "../src/headers/clusterMix.h"
#include "gtest/gtest.h"

class ClusterMixTest : public ::testing::Test {
protected:
    // std::string filename;
    // Coal::EventsMap all_events{};
    // Coal::EventsMap Cell{};
    // const Coal::ClusterParams he4_params = {
    //         4,
    //         20,
    //         4,
    //         {1.45, 1.45, 1.45},
    //         1. / 16,
    //         444,
    //         {2212, 2212, 2112, 2112},
    //         1.0e-8,
    //         {{-0.8, 0.8}, {-0.8, 0.8}, {-0.8, 0.8}, {-0.8, 0.8}},
    //         {{-0.5, 0.5}},
    //         {0.2, 25}};
    //
    // const Coal::ClusterParams deutron_params = {2,
    //                                             100,
    //                                             10,
    //                                             {1.96},
    //                                             3. / 4,
    //                                             222,
    //                                             {2212, 2112},
    //                                             0.00001,
    //                                             {{-0.8, 0.8}, {-0.8, 0.8}},
    //                                             {{-0.5, 0.5}},
    //                                             {0.2, 20}};
    //
    // std::vector<Coal::ParticleArray> subCell{};
    //
    // void SetUp() override {
    //     filename   = "../../../input/100/particle_lists.oscar";
    //     all_events = Coal::Node::readFile(filename, "smash");
    //
    //     Coal::RandomNumber::getInstance().seed(123456);
    //     Cell    = Coal::selectEvents(all_events, he4_params);
    //     subCell = Coal::selectParticles(Cell, he4_params);
    // }
};

TEST_F(ClusterMixTest, indexToMultiIndexWorks) {

    const std::vector counts = {152, 147, 160};

    //Test indexToMultiIndex
    const std::vector<int> result1 = Coal::indexToMultiIndex(23840, counts);
    const std::vector expected1    = {0, 0, 1, 0};
    ASSERT_EQ(expected1, result1);
    for (auto i = 0; i < counts.size(); ++i) {
        ASSERT_LE(result1[i], counts[i]);
    }
}
//
TEST_F(ClusterMixTest, multiIndexToIndexWorks) {
    const std::vector counts = {152, 147, 160};
    const std::vector index  = {1, 2, 0};
    const long long result   = Coal::multiIndexToIndex(index, counts);
    const long long expected = 174;
    ASSERT_EQ(expected, result);
    std::cout << "result:" << result << "except:" << expected << std::endl;
}

TEST_F(ClusterMixTest, jumpValidLoopWorks) {
    const std::vector counts  = {304, 279, 373, 370};
    const std::vector result1 = {300, 278, 330, 300};

    const std::vector result = Coal::jumpValidLoop(result1, counts, 1);

    const std::vector expected{303, 278, 372, 360};
    ASSERT_EQ(expected, result);
    for (auto i = 0; i < counts.size(); ++i) {
        std::cout << result[i] << ",";
    }
}


// TEST_F(ClusterMixTest, Nbody) {
//     const int N = he4_params.NBody;
//
//     auto counts       = std::vector<int>(N, 0);
//     long long loopMax = 1;
//     for (auto i = 0; i < N; ++i) {
//         counts[i] = static_cast<int>(subCell[i].size());
//         std::cout << counts[i] << ",";
//         loopMax *= counts[i];
//     }
//     std::cout << "loop:" << loopMax << std::endl;
//     // const std::vector muliIndex = {1, 0, 17, 42};
//     const std::vector<int> muliIndex = Coal::indexToMultiIndex(4653498, counts);
//     for (auto i = 0; i < N; ++i) {
//         std::cout << muliIndex[i] << ",";
//     }
//     Coal::ParticleArray particles{};
//     for (auto i = 0; i < N; ++i) {
//         particles.push_back(subCell[i][muliIndex[i]]);
//     }
//
//     auto [disx, dispp] = Coal::JacobiCoordinates(particles);
//     for (auto i = 0; i < N - 1; ++i) {
//         std::cout << "dx:" << disx[i] << " "
//                   << "dp:" << dispp[i] << std::endl;
//     }
//     const auto result   = loopMax;
//     const auto expected = 1;
//     ASSERT_EQ(expected, result);
// }
