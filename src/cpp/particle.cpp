//
// Created by mafu on 1/15/2024.
//
#include "../headers/particle.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <ranges>


bool Coal::Particle::operator!=(const Particle &other) const {
    constexpr double epsilon = 1e-6;
    return std::abs(t - other.t) > epsilon || std::abs(x - other.x) > epsilon ||
           std::abs(y - other.y) > epsilon || std::abs(z - other.z) > epsilon ||
           std::abs(p0 - other.p0) > epsilon ||
           std::abs(px - other.px) > epsilon ||
           std::abs(py - other.py) > epsilon ||
           std::abs(pz - other.pz) > epsilon;
}

std::tuple<double, double, double>
Coal::Particle::get4Momentum(const std::vector<Particle> &particles) {
    return std::accumulate(
            particles.begin(), particles.end(), std::tuple{0.0, 0.0, 0.0},
            [](const auto &a, const auto &b) {
                auto [px, py, pz] = a;
                return std::tuple{px + b.px, py + b.py, pz + b.pz};
            });
}

std::tuple<double, double, double, double>
Coal::Particle::getPosition(const std::vector<Particle> &particles) {
    return std::accumulate(
            particles.begin(), particles.end(), std::tuple{0.0, 0.0, 0.0, 0.0},
            [](const auto &a, const auto &b) {
                auto [x1, y1, z1, m1] = a;
                return std::tuple{x1 + b.x, y1 + b.y, z1 + b.z, m1 + b.mass};
            });
}


void Coal::Particle::getResultParticleData(
        const std::vector<Particle> &particles) {

    const int n = static_cast<int>(particles.size());
    if (n == 0)
        return;
    charge = 0;
    for (auto &particle: particles) {
        charge += particle.charge;
    }

    freeze_out_time = std::ranges::max_element(particles, [](const auto &a,
                                                             const auto &b) {
                          return a.freeze_out_time < b.freeze_out_time;
                      })->freeze_out_time;

    auto [sum_x, sum_y, sum_z, sum_mass] = getPosition(particles);
    auto [sum_px, sum_py, sum_pz]        = get4Momentum(particles);
    t                                    = freeze_out_time;
    x                                    = sum_x / n;
    y                                    = sum_y / n;
    z                                    = sum_z / n;
    px                                   = sum_px;
    py                                   = sum_py;
    pz                                   = sum_pz;
    mass                                 = sum_mass;
    // p0 = std::sqrt(px * px + py * py + pz * pz + mass * mass);
    p0 = std::accumulate(particles.begin(), particles.end(), 0.0,
                         [](const auto &a, const auto &b) { return a + b.p0; });
}

double Coal::Particle::getRapidity() const {
    const double E = std::sqrt(px * px + py * py + pz * pz + mass * mass);
    if (E <= pz) {
        return 0.0;
    }
    return 0.5 * std::log((E + pz) / (E - pz));
}

double Coal::Particle::getPT() const { return std::sqrt(px * px + py * py); }

double Coal::Particle::getArtifactRapidity() const {
    const double p_total = std::sqrt(px * px + py * py + pz * pz);
    if (p_total <= pz) {
        return 0.0;
    }
    return 0.5 * std::log((p_total + pz) / (p_total - pz));
}

void Coal::Particle::getFreezeOutPosition() {
    if (freeze_out_time == 0.0 || std::abs(t - freeze_out_time) < 1e-6 ||
        t == 0.0) {
    } else {
        const double vx = px / p0;
        const double vy = py / p0;
        const double vz = pz / p0;
        const double dt = t - freeze_out_time;
        x -= dt * vx;
        y -= dt * vy;
        z -= dt * vz;
    }
}

void Coal::Particle::updatePosition(const double t_max) {
    const double dt = t_max - freeze_out_time;
    x += dt * px / p0;
    y += dt * py / p0;
    z += dt * pz / p0;
}

Coal::Particle Coal::Particle::lorentzBoost(const double beta_x,
                                            const double beta_y,
                                            const double beta_z) const {
    Particle boost_p = *this;
    double beta2     = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;
    if (beta2 == 0.0 || beta2 < 1e-6) {
        return boost_p;
    }
    if (beta2 > 0.9999999) {
        beta2 = 0.9999999;
    }
    const double gamma   = 1.0 / sqrt(1.0 - beta2);
    const double x_lam00 = gamma;
    const double x_lam01 = -gamma * beta_x;
    const double x_lam02 = -gamma * beta_y;
    const double x_lam03 = -gamma * beta_z;
    const double x_lam11 = 1.0 + (gamma - 1.0) * beta_x * beta_x / beta2;
    const double x_lam22 = 1.0 + (gamma - 1.0) * beta_y * beta_y / beta2;
    const double x_lam33 = 1.0 + (gamma - 1.0) * beta_z * beta_z / beta2;
    const double x_lam12 = (gamma - 1.0) * beta_x * beta_y / beta2;
    const double x_lam13 = (gamma - 1.0) * beta_x * beta_z / beta2;
    const double x_lam23 = (gamma - 1.0) * beta_y * beta_z / beta2;
    const double new_t =
            freeze_out_time * x_lam00 + x * x_lam01 + y * x_lam02 + z * x_lam03;
    const double new_x =
            freeze_out_time * x_lam01 + x * x_lam11 + y * x_lam12 + z * x_lam13;
    const double new_y =
            freeze_out_time * x_lam02 + x * x_lam12 + y * x_lam22 + z * x_lam23;
    const double new_z =
            freeze_out_time * x_lam03 + x * x_lam13 + y * x_lam23 + z * x_lam33;
    const double new_p0 =
            p0 * x_lam00 + px * x_lam01 + py * x_lam02 + pz * x_lam03;
    const double new_px =
            p0 * x_lam01 + px * x_lam11 + py * x_lam12 + pz * x_lam13;
    const double new_py =
            p0 * x_lam02 + px * x_lam12 + py * x_lam22 + pz * x_lam23;
    const double new_pz =
            p0 * x_lam03 + px * x_lam13 + py * x_lam23 + pz * x_lam33;
    boost_p.freeze_out_time = new_t;
    boost_p.x               = new_x;
    boost_p.y               = new_y;
    boost_p.z               = new_z;
    boost_p.p0              = new_p0;
    boost_p.px              = new_px;
    boost_p.py              = new_py;
    boost_p.pz              = new_pz;
    return boost_p;
}
