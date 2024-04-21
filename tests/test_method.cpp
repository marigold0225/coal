//
// Created by mafu on 1/6/2024.
//
#include "method.h"
#include "random.h"
#include "flow.h"
#include "shortcut.h"
#include "gtest/gtest.h"
#include <bits/fs_fwd.h>
#include <bits/fs_path.h>
#include <iostream>

class MethodTest : public testing::Test {
protected:
    std::string filename;
    Coal::EventsMap allEvents;
    Coal::EventsMap Cell;
    std::vector<Eigen::MatrixXd> subCell;
    Coal::ParticleEventMap proton;
    Coal::ParticleEventMap neutron;
    Coal::ParticleArray particle_before{};

    void SetUp() override {
        filename = "../../../coal/input/100/particle_lists.oscar";
        // Assuming load returns two vectors of Coal::Particle
        std::cout << "currentPath is:" << std::filesystem::current_path()
                  << std::endl;
        YAML::Node config  = YAML::LoadFile("../../input/input.yaml");
        auto ClusterParams = config["ClusterParams"];
        const Coal::ClusterParams cluster(ClusterParams["Helium4"]);


        Coal::RandomNumber::getInstance().seed(123456);
        allEvents = Coal::FileLoader::readFile(filename, "smash");
        auto resolution = Coal::getResolutionMap(allEvents, 1);

        Cell    = selectEvents(allEvents, cluster, resolution, 1, TODO);
        subCell = selectParticles(Cell, cluster);
    }
};

TEST_F(MethodTest, JacobiFunctionWorks) {
    const std::vector massArray_test1 = {1.0, 2.0};
    const std::vector massArray_test2 = {2.0, 3.0, 5.0};
    const std::vector massArray_test3 = {1.0, 1.0, 2.0, 3.0};
    const Eigen::MatrixXd result1     = Coal::MassMatrix(massArray_test1);
    const Eigen::MatrixXd result2     = Coal::MassMatrix(massArray_test2);
    const Eigen::MatrixXd result3     = Coal::MassMatrix(massArray_test3);

    Eigen::MatrixXd expected1(2, 2);
    expected1 << 0.3333, 0.6667, 0.7071, -0.7071;

    Eigen::MatrixXd expected2(3, 3);
    expected2 << 0.2000, 0.3000, 0.5000, 0.7071, -0.7071, 0, 0.3266, 0.4899,
            -0.8165;

    Eigen::MatrixXd expected3(4, 4);
    expected3 << 0.1429, 0.1429, 0.2857, 0.4286, 0.7071, -0.7071, 0, 0, 0.4082,
            0.4082, -0.8165, 0, 0.2165, 0.2165, 0.4330, -0.8660;

    if (!expected1.isApprox(result1, 1e-4)) {
        std::cout << "Result 1: " << std::endl << result1 << std::endl;
        std::cout << "Expected 1: " << std::endl << expected1 << std::endl;
    }
    ASSERT_TRUE(expected1.isApprox(result1, 1e-4));

    if (!expected2.isApprox(result2, 1e-4)) {
        std::cout << "Result 2: " << std::endl << result2 << std::endl;
        std::cout << "Expected 2: " << std::endl << expected2 << std::endl;
    }
    ASSERT_TRUE(expected2.isApprox(result2, 1e-4));

    if (!expected3.isApprox(result3, 1e-4)) {
        std::cout << "Result 3: " << std::endl << result3 << std::endl;
        std::cout << "Expected 3: " << std::endl << expected3 << std::endl;
    }
    ASSERT_TRUE(expected3.isApprox(result3, 1e-4));
}

// TEST_F(MethodTest, lorentzboost) {
//     Coal::Particle result{};
//     Coal::Particle result2{};
//     const Coal::Particle proton1  = subCell[0][0];
//     const Coal::Particle proton2  = subCell[1][0];
//     const Coal::Particle neutron1 = subCell[2][0];
//     const Coal::Particle neutron2 = subCell[3][0];
//     const std::vector particles   = {proton1, proton2, neutron1, neutron2};
//
//     result.getResultParticleData(particles);
//     std::cout << "x:" << result.x << "y:" << result.y << "z:" << result.z
//               << "t:" << result.freeze_out_time << "\n"
//               << "px:" << result.px << "py:" << result.py << "pz:" << result.pz
//               << "p0:" << result.p0 << "\n";
//     //
//
//     //
//     double total_px = 0.0, total_py = 0.0, total_pz = 0.0, total_p0 = 0.0;
//
//     for (auto &particle: particles) {
//         total_px += particle.px;
//         total_py += particle.py;
//         total_pz += particle.pz;
//         total_p0 += particle.p0;
//     }
//
//
//     const double beta_x = total_px / total_p0;
//     const double beta_y = total_py / total_p0;
//     const double beta_z = total_pz / total_p0;
//
//
//     Coal::ParticleArray boostparticles;
//     for (auto &particle: particles) {
//         boostparticles.push_back(particle.lorentzBoost(beta_x, beta_y, beta_z));
//     }
//     const double t_max =
//             std::ranges::max_element(boostparticles, [](const auto &a,
//                                                         const auto &b) {
//                 return a.freeze_out_time < b.freeze_out_time;
//             })->freeze_out_time;
//     for (auto &particle: boostparticles) {
//         particle.updatePosition(t_max);
//     }
//     result2.getResultParticleData(boostparticles);
//     const auto result3 = result2.lorentzBoost(-beta_x, -beta_y, -beta_z);
//     std::cout << "x:" << result3.x << "y:" << result3.y << "z:" << result3.z
//               << "t:" << result3.freeze_out_time << "\n"
//               << "px:" << result3.px << "py:" << result3.py
//               << "pz:" << result3.pz << "p0:" << result3.p0 << "\n";
//     ASSERT_NEAR(result.x, result3.px, 1.0);
// }