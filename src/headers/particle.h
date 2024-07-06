//
// Created by mafu on 1/15/2024.
//
#pragma once
// #define EIGEN_USE_MKL_ALL
#include <vector>
#include <Eigen/Dense>
#define M_PI 3.14159265358979323846

namespace Coal {
    struct Particle {
        double t, x, y, z;
        double p0, px, py, pz;
        double mass;
        int pdg;
        double freeze_out_time;
        [[maybe_unused]] int charge;
        [[maybe_unused]] double probability;

        // math functions
        bool operator!=(const Particle &other) const;

        static std::tuple<double, double, double>
        get4Momentum(const std::vector<Particle> &particles);

        static std::tuple<double, double, double, double>
        getPosition(const std::vector<Particle> &particles);

        void getResultParticleData(const std::vector<Particle> &particles);

        void updatePosition(double t_max);

        [[nodiscard]] Particle lorentzBoost(double beta_x, double beta_y,
                                            double beta_z) const;

        [[nodiscard]] double getRapidity() const;

        [[nodiscard]] double getPT() const;

        [[nodiscard]] double getPhi() const;

        [[nodiscard]] double getPseudoRapidity() const;

        void getFreezeOutPosition();

        [[nodiscard]] Eigen::MatrixXd convertToMatrix() const;
    };
    std::ostream& operator<<(std::ostream& os, const Particle& particle);
}// namespace Coal
