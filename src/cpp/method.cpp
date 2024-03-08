//
// Created by mafu on 1/8/2024.
//
#include "../headers/method.h"
#include "omp.h"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <ranges>

Coal::ClusterParams::ClusterParams(const YAML::Node &node) {
    NBody       = node["NBody"].as<int>();
    Loop        = node["Loop"].as<int>();
    mixEvents   = node["MixEvents"].as<int>();
    nmix        = mixEvents / NBody;
    eventFactor = Loop * pow(nmix, NBody);
    probFactor  = pow(8, NBody - 1);
    precision   = node["Precision"].as<int>();
    if (const auto gc_value = node["gc"].as<std::string>();
        gc_value.find('/') != std::string::npos) {
        gc = parseFraction(gc_value);
    } else {
        gc = node["gc"].as<double>();
    }
    pdg            = node["PDG"].as<int>();
    probabilityCut = node["ProbCut"].as<double>();
    MassArray      = node["MassArray"].as<std::vector<double>>();
    SigArray       = node["Sig"].as<std::vector<double>>();
    PDGArray       = node["From"].as<std::vector<int>>();
    ptBins         = node["Pt"].as<std::pair<double, int>>();
    originRapidity = node["RapidityCut"].as<RapidityArray>();
    targetRapidity = node["TargetRapidity"].as<RapidityRange>();
    setMassArray(MassArray);
}
void Coal::ClusterParams::setMassArray(const std::vector<double> &massArray) {
    M.resize(massArray.size() - 1);
    M_inv_t.resize(massArray.size() - 1);
    for (int i = 1; i < massArray.size(); ++i) {
        std::vector temMassArray(massArray.begin(), massArray.begin() + i + 1);
        M[i - 1]       = MassMatrix(temMassArray);
        M_inv_t[i - 1] = M[i - 1].inverse().transpose();
    }
}

double Coal::parseFraction(const std::string &fraction) {
    std::istringstream iss(fraction);
    double numerator, denominator;
    if (char slash; !(iss >> numerator >> slash >> denominator)) {
        std::cerr << "Error: cannot parse fraction " << fraction << std::endl;
        exit(1);
    }
    return numerator / denominator;
}
Coal::RapidityArray Coal::defineRapidityRange(const YAML::Node &node) {
    RapidityArray rapidityArray;
    for (const auto rapidity =
                 node["Output"]["RapidityRange"]
                         .as<std::vector<std::pair<double, double>>>();
         const auto &rap: rapidity) {
        rapidityArray.emplace_back(rap);
    }
    return rapidityArray;
}
bool Coal::checkFileExits(const std::string &filename,
                          const std::vector<std::string> &labels,
                          const std::string &fileType) {
    for (auto &label: labels) {
        std::string fileName = filename;
        fileName.append("/")
                .append(label)
                .append("/")
                .append(fileType)
                .append("_")
                .append(label)
                .append(".dat");
        if (!std::filesystem::exists(fileName)) {
            return false;
        }
    }
    return true;
}
std::string Coal::constructFilename(const std::string &filename,
                                    const std::string &fileType,
                                    const std::string &label) {
    return filename + "/" + label + "/" + fileType + "_" + label + ".dat";
}
std::vector<double> Coal::linspace(const double start, const double end,
                                   const int num) {
    std::vector<double> Linspaced;
    const double delta = (end - start) / (num - 1);

    Linspaced.reserve(num - 1);
    for (int i = 0; i < num - 1; ++i) {
        Linspaced.push_back(start + delta * i);
    }
    Linspaced.push_back(end);

    return Linspaced;
}

Eigen::MatrixXd Coal::MassMatrix(const std::vector<double> &massArray) {
    const int n       = static_cast<int>(massArray.size());
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);

    const double totalMass = Eigen::VectorXd::Map(massArray.data(), n).sum();

    M.row(0) = Eigen::VectorXd::Map(massArray.data(), n) / totalMass;

    for (int k = 1; k < n; ++k) {
        if (k == 1) {
            M(k, 0) = 1.0 / sqrt(2);
            M(k, 1) = -1.0 / sqrt(2);
        } else {
            double sum = 0;
            for (int i = 0; i < k; ++i) {
                sum += massArray[i];
            }
            for (int i = 0; i < k; ++i) {
                M(k, i) = sqrt(static_cast<double>(k) / (k + 1)) *
                          massArray[i] / sum;
            }
            M(k, k) = -sqrt(static_cast<double>(k) / (k + 1));
        }
    }

    return M;
}

Coal::ParticleArray Coal::boostToCOM(const ParticleArray &particles,
                                     const Particle &targetParticle) {
    using std::ranges::max_element;
    const double beta_x = targetParticle.px / targetParticle.p0;
    const double beta_y = targetParticle.py / targetParticle.p0;
    const double beta_z = targetParticle.pz / targetParticle.p0;

    ParticleArray boostedParticles;
    for (auto &particle: particles) {
        boostedParticles.push_back(
                particle.lorentzBoost(beta_x, beta_y, beta_z));
    }

    const double t_max =
            max_element(boostedParticles, [](const auto &a, const auto &b) {
                return a.freeze_out_time < b.freeze_out_time;
            })->freeze_out_time;

    for (auto &particle: boostedParticles) {
        particle.updatePosition(t_max);
    }
    return boostedParticles;
}
void Coal::lorenzBoostMatrix(Eigen::MatrixXd &particles, const double beta_x,
                             const double beta_y, const double beta_z) {
    double beta2       = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;
    const double gamma = 1.0 / sqrt(1.0 - beta2);
    if (beta2 == 0 || beta2 < 1e-6) {
        return;
    }
    if (beta2 > 0.99999) {
        beta2 = 0.99999;
    }
    Eigen::Matrix4d lambda;
    lambda << gamma, -gamma * beta_x, -gamma * beta_y, -gamma * beta_z,
            -gamma * beta_x, 1 + (gamma - 1) * beta_x * beta_x / beta2,
            (gamma - 1) * beta_x * beta_y / beta2,
            (gamma - 1) * beta_x * beta_z / beta2, -gamma * beta_y,
            (gamma - 1) * beta_y * beta_x / beta2,
            1 + (gamma - 1) * beta_y * beta_y / beta2,
            (gamma - 1) * beta_y * beta_z / beta2, -gamma * beta_z,
            (gamma - 1) * beta_z * beta_x / beta2,
            (gamma - 1) * beta_z * beta_y / beta2,
            1 + (gamma - 1) * beta_z * beta_z / beta2;

    for (int i = 0; i < particles.rows(); ++i) {
        Eigen::Vector4d particleMomentum(particles(i, 5), particles(i, 1),
                                         particles(i, 2), particles(i, 3));
        Eigen::Vector4d particlePosition(particles(i, 9), particles(i, 6),
                                         particles(i, 7), particles(i, 8));

        Eigen::Vector4d boostedMomentum = lambda * particleMomentum;
        Eigen::Vector4d boostedPosition = lambda * particlePosition;

        particles(i, 1) = boostedMomentum(1);
        particles(i, 2) = boostedMomentum(2);
        particles(i, 3) = boostedMomentum(3);
        particles(i, 5) = boostedMomentum(0);
        particles(i, 6) = boostedPosition(1);
        particles(i, 7) = boostedPosition(2);
        particles(i, 8) = boostedPosition(3);
        particles(i, 9) = boostedPosition(0);
    }
}
Eigen::MatrixXd Coal::boostToComMatrix(const Eigen::MatrixXd &particles) {
    const double beta_x = particles.col(1).sum() / particles.col(5).sum();
    const double beta_y = particles.col(2).sum() / particles.col(5).sum();
    const double beta_z = particles.col(3).sum() / particles.col(5).sum();

    Eigen::MatrixXd boostedParticles = particles;

    lorenzBoostMatrix(boostedParticles, beta_x, beta_y, beta_z);

    const double t_max = boostedParticles.col(9).maxCoeff();
    for (int i = 0; i < boostedParticles.rows(); ++i) {
        boostedParticles(i, 6) += (t_max - boostedParticles(i, 9)) *
                                  boostedParticles(i, 1) /
                                  boostedParticles(i, 5);
        boostedParticles(i, 7) += (t_max - boostedParticles(i, 9)) *
                                  boostedParticles(i, 2) /
                                  boostedParticles(i, 5);
        boostedParticles(i, 8) += (t_max - boostedParticles(i, 9)) *
                                  boostedParticles(i, 3) /
                                  boostedParticles(i, 5);
    }
    return boostedParticles;
}
std::vector<Eigen::Matrix4d>
Coal::calculateLorentz(const Eigen::VectorXd &betaX,
                       const Eigen::VectorXd &betaY,
                       const Eigen::VectorXd &betaZ) {
    std::vector<Eigen::Matrix4d> lorentzMatrices;
    lorentzMatrices.reserve(betaX.size());
    for (int i = 0; i < betaX.size(); ++i) {
        const double beta_x = betaX(i);
        const double beta_y = betaY(i);
        const double beta_z = betaZ(i);
        double beta2 = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;

        if (beta2 == 0 || beta2 < 1e-6) {
            lorentzMatrices.emplace_back(Eigen::Matrix4d::Identity());
            continue;
        }
        if (beta2 > 0.99999) {
            beta2 = 0.99999;
        }

        double gamma = 1.0 / sqrt(1.0 - beta2);

        Eigen::Matrix4d lambda;
        lambda << gamma, -gamma * beta_x, -gamma * beta_y, -gamma * beta_z,
                -gamma * beta_x, 1 + (gamma - 1) * beta_x * beta_x / beta2,
                (gamma - 1) * beta_x * beta_y / beta2,
                (gamma - 1) * beta_x * beta_z / beta2, -gamma * beta_y,
                (gamma - 1) * beta_y * beta_x / beta2,
                1 + (gamma - 1) * beta_y * beta_y / beta2,
                (gamma - 1) * beta_y * beta_z / beta2, -gamma * beta_z,
                (gamma - 1) * beta_z * beta_x / beta2,
                (gamma - 1) * beta_z * beta_y / beta2,
                1 + (gamma - 1) * beta_z * beta_z / beta2;

        lorentzMatrices.push_back(lambda);
    }
    return lorentzMatrices;
}

void Coal::applyLorentzBoost(Eigen::Ref<Eigen::MatrixXd> combinedX,
                             Eigen::Ref<Eigen::MatrixXd> combinedY,
                             Eigen::Ref<Eigen::MatrixXd> combinedZ,
                             Eigen::Ref<Eigen::MatrixXd> combinedT,
                             Eigen::Ref<Eigen::MatrixXd> combinedPX,
                             Eigen::Ref<Eigen::MatrixXd> combinedPY,
                             Eigen::Ref<Eigen::MatrixXd> combinedPZ,
                             Eigen::Ref<Eigen::MatrixXd> combinedP0,
                             const std::vector<Eigen::Matrix4d> &lorentz) {

    Eigen::Vector4d particleMomentum;
    Eigen::Vector4d particlePosition;
    Eigen::Vector4d boostedMomentum;
    Eigen::Vector4d boostedPosition;

    for (int i = 0; i < combinedX.rows(); ++i) {
        const Eigen::Matrix4d &lorentzMatrix = lorentz[i];

        for (int j = 0; j < combinedX.cols(); ++j) {

            particleMomentum =
                    Eigen::Vector4d(combinedP0(i, j), combinedPX(i, j),
                                    combinedPY(i, j), combinedPZ(i, j));
            particlePosition =
                    Eigen::Vector4d(combinedT(i, j), combinedX(i, j),
                                    combinedY(i, j), combinedZ(i, j));

            boostedMomentum.noalias() = lorentzMatrix * particleMomentum;
            boostedPosition.noalias() = lorentzMatrix * particlePosition;

            combinedPX(i, j) = boostedMomentum(1);
            combinedPY(i, j) = boostedMomentum(2);
            combinedPZ(i, j) = boostedMomentum(3);
            combinedP0(i, j) = boostedMomentum(0);
            combinedX(i, j)  = boostedPosition(1);
            combinedY(i, j)  = boostedPosition(2);
            combinedZ(i, j)  = boostedPosition(3);
            combinedT(i, j)  = boostedPosition(0);
        }
    }
}
void Coal::applyLorentzBoostCols(Eigen::Ref<Eigen::MatrixXd> combinedX,
                                 Eigen::Ref<Eigen::MatrixXd> combinedY,
                                 Eigen::Ref<Eigen::MatrixXd> combinedZ,
                                 Eigen::Ref<Eigen::MatrixXd> combinedT,
                                 Eigen::Ref<Eigen::MatrixXd> combinedPX,
                                 Eigen::Ref<Eigen::MatrixXd> combinedPY,
                                 Eigen::Ref<Eigen::MatrixXd> combinedPZ,
                                 Eigen::Ref<Eigen::MatrixXd> combinedP0,
                                 const std::vector<Eigen::Matrix4d> &lorentz) {

    for (int i = 0; i < combinedX.rows(); ++i) {
        const Eigen::Matrix4d &lorentzMatrix = lorentz[i];

        for (int j = 0; j < combinedX.cols(); ++j) {
            // 提取矩阵和向量的元素
            double p0 = combinedP0(i, j);
            double px = combinedPX(i, j);
            double py = combinedPY(i, j);
            double pz = combinedPZ(i, j);
            double t  = combinedT(i, j);
            double x  = combinedX(i, j);
            double y  = combinedY(i, j);
            double z  = combinedZ(i, j);

            // 展开乘法计算
            combinedP0(i, j) =
                    lorentzMatrix(0, 0) * p0 + lorentzMatrix(0, 1) * px +
                    lorentzMatrix(0, 2) * py + lorentzMatrix(0, 3) * pz;
            combinedPX(i, j) =
                    lorentzMatrix(1, 0) * p0 + lorentzMatrix(1, 1) * px +
                    lorentzMatrix(1, 2) * py + lorentzMatrix(1, 3) * pz;
            combinedPY(i, j) =
                    lorentzMatrix(2, 0) * p0 + lorentzMatrix(2, 1) * px +
                    lorentzMatrix(2, 2) * py + lorentzMatrix(2, 3) * pz;
            combinedPZ(i, j) =
                    lorentzMatrix(3, 0) * p0 + lorentzMatrix(3, 1) * px +
                    lorentzMatrix(3, 2) * py + lorentzMatrix(3, 3) * pz;

            combinedT(i, j) = lorentzMatrix(0, 0) * t +
                              lorentzMatrix(0, 1) * x +
                              lorentzMatrix(0, 2) * y + lorentzMatrix(0, 3) * z;
            combinedX(i, j) = lorentzMatrix(1, 0) * t +
                              lorentzMatrix(1, 1) * x +
                              lorentzMatrix(1, 2) * y + lorentzMatrix(1, 3) * z;
            combinedY(i, j) = lorentzMatrix(2, 0) * t +
                              lorentzMatrix(2, 1) * x +
                              lorentzMatrix(2, 2) * y + lorentzMatrix(2, 3) * z;
            combinedZ(i, j) = lorentzMatrix(3, 0) * t +
                              lorentzMatrix(3, 1) * x +
                              lorentzMatrix(3, 2) * y + lorentzMatrix(3, 3) * z;
        }
    }
}


std::tuple<std::vector<double>, std::vector<double>>
Coal::JacobiCoordinatesMatrix(const Eigen::MatrixXd &particles,
                              const ClusterParams &params) {
    const auto N = particles.rows();
    Eigen::VectorXd x(N), y(N), z(N), px(N), py(N), pz(N);

    x  = particles.col(6);
    y  = particles.col(7);
    z  = particles.col(8);
    px = particles.col(1);
    py = particles.col(2);
    pz = particles.col(3);

    auto transformCoordinates = [&](const Eigen::VectorXd &vec,
                                    const Eigen::MatrixXd &mat) {
        return (mat * vec).tail(vec.size() - 1);
    };

    Eigen::VectorXd x_jacobi  = transformCoordinates(x, params.M[N - 2]);
    Eigen::VectorXd y_jacobi  = transformCoordinates(y, params.M[N - 2]);
    Eigen::VectorXd z_jacobi  = transformCoordinates(z, params.M[N - 2]);
    Eigen::VectorXd px_jacobi = transformCoordinates(px, params.M_inv_t[N - 2]);
    Eigen::VectorXd py_jacobi = transformCoordinates(py, params.M_inv_t[N - 2]);
    Eigen::VectorXd pz_jacobi = transformCoordinates(pz, params.M_inv_t[N - 2]);

    std::vector<double> d_r(N - 1), d_p(N - 1);
    for (auto i = 0; i < N - 1; ++i) {
        d_r[i] = sqrt(x_jacobi(i) * x_jacobi(i) + y_jacobi(i) * y_jacobi(i) +
                      z_jacobi(i) * z_jacobi(i));
        d_p[i] =
                sqrt(px_jacobi(i) * px_jacobi(i) + py_jacobi(i) * py_jacobi(i) +
                     pz_jacobi(i) * pz_jacobi(i));
    }
    return std::make_tuple(d_r, d_p);
}

std::tuple<std::vector<double>, std::vector<double>>
Coal::JacobiCoordinates(const ParticleArray &particles,
                        const ClusterParams &params) {
    const auto N = static_cast<int>(particles.size());

    Eigen::VectorXd x(N), y(N), z(N), px(N), py(N), pz(N);

    for (auto i = 0; i < N; ++i) {
        x(i)  = particles[i].x;
        y(i)  = particles[i].y;
        z(i)  = particles[i].z;
        px(i) = particles[i].px;
        py(i) = particles[i].py;
        pz(i) = particles[i].pz;
    }

    auto transformCoordinates = [&](const Eigen::VectorXd &vec,
                                    const Eigen::MatrixXd &mat) {
        return (mat * vec).tail(vec.size() - 1);
    };

    Eigen::VectorXd x_jacobi  = transformCoordinates(x, params.M[N - 2]);
    Eigen::VectorXd y_jacobi  = transformCoordinates(y, params.M[N - 2]);
    Eigen::VectorXd z_jacobi  = transformCoordinates(z, params.M[N - 2]);
    Eigen::VectorXd px_jacobi = transformCoordinates(px, params.M_inv_t[N - 2]);
    Eigen::VectorXd py_jacobi = transformCoordinates(py, params.M_inv_t[N - 2]);
    Eigen::VectorXd pz_jacobi = transformCoordinates(pz, params.M_inv_t[N - 2]);

    std::vector<double> d_r(N - 1), d_p(N - 1);
    for (auto i = 0; i < N - 1; ++i) {
        d_r[i] = sqrt(x_jacobi(i) * x_jacobi(i) + y_jacobi(i) * y_jacobi(i) +
                      z_jacobi(i) * z_jacobi(i));
        d_p[i] =
                sqrt(px_jacobi(i) * px_jacobi(i) + py_jacobi(i) * py_jacobi(i) +
                     pz_jacobi(i) * pz_jacobi(i));
    }
    return std::make_tuple(d_r, d_p);
}

Coal::EventsMap Coal::selectEvents(const EventsMap &eventMap,
                                   const ClusterParams &params) {
    EventsMap sampleEvents;

    int selectEventsSize = params.mixEvents;

    auto &generator = RandomNumber::getInstance().getGenerator();

    std::vector<int> eventIDList;

    for (const auto &key: eventMap | std::views::keys) {
        eventIDList.push_back(key);
    }

    selectEventsSize =
            std::min(selectEventsSize, static_cast<int>(eventIDList.size()));

    std::ranges::shuffle(eventIDList, generator);
    for (int i = 0; i < selectEventsSize; ++i) {
        auto eventID          = eventIDList[i];
        sampleEvents[eventID] = eventMap.at(eventID);
    }
    return sampleEvents;
}
Coal::MultiParticleArray Coal::selectParticles(const EventsMap &eventMap,
                                               const ClusterParams &params) {

    const auto eachPDGEvent = params.nmix;
    const auto PDGs         = params.PDGArray;
    MultiParticleArray result;
    std::vector<int> eventIDList;
    for (const auto &ID: eventMap | std::views::keys) {
        eventIDList.push_back(ID);
    }

    auto &generator = RandomNumber::getInstance().getGenerator();
    std::ranges::shuffle(eventIDList, generator);
    // std::ranges::reverse(eventIDList);
    ParticleArray particlesForThisPDG{};

    for (size_t i = 0; i < PDGs.size(); ++i) {
        int pdgCode = PDGs[i];
        int count   = 0;

        while (count < eachPDGEvent && !eventIDList.empty()) {
            int eventID = eventIDList.back();
            eventIDList.pop_back();

            if (auto eventIt = eventMap.find(eventID);
                eventIt != eventMap.end()) {
                if (auto particleIt = eventIt->second.find(pdgCode);
                    particleIt != eventIt->second.end()) {
                    for (const auto &particle: particleIt->second) {
                        if (const double rapidity = particle.getRapidity();
                            rapidity > params.originRapidity[i].first &&
                            rapidity < params.originRapidity[i].second) {
                            particlesForThisPDG.push_back(particle);
                        }
                    }
                } else {
                    std::cout << "pdgCode: " << pdgCode
                              << " not found in eventID: " << eventID
                              << std::endl;
                    return result;
                }
            }
            ++count;
        }
        std::ranges::shuffle(particlesForThisPDG, generator);
        result.push_back(std::move(particlesForThisPDG));
        particlesForThisPDG.clear();
    }
    return result;
}
std::vector<Eigen::MatrixXd>
Coal::selectParticlesMatrix(const EventsMap &eventMap,
                            const ClusterParams &params) {
    const auto eachPDGEvent = params.nmix;
    const auto PDGs         = params.PDGArray;
    std::vector<Eigen::MatrixXd> result;
    std::vector<int> eventIDList;
    for (const auto &ID: eventMap | std::views::keys) {
        eventIDList.push_back(ID);
    }

    auto &generator = RandomNumber::getInstance().getGenerator();
    std::ranges::shuffle(eventIDList, generator);
    // std::ranges::reverse(eventIDList);
    ParticleArray particlesForThisPDG{};
    Eigen::MatrixXd matrix(0, 11);
    for (size_t i = 0; i < PDGs.size(); ++i) {
        int pdgCode = PDGs[i];
        int count   = 0;

        while (count < eachPDGEvent && !eventIDList.empty()) {
            int eventID = eventIDList.back();
            eventIDList.pop_back();

            if (auto eventIt = eventMap.find(eventID);
                eventIt != eventMap.end()) {
                if (auto particleIt = eventIt->second.find(pdgCode);
                    particleIt != eventIt->second.end()) {
                    for (const auto &particle: particleIt->second) {
                        if (const double rapidity = particle.getRapidity();
                            rapidity > params.originRapidity[i].first &&
                            rapidity < params.originRapidity[i].second) {
                            particlesForThisPDG.push_back(particle);
                        }
                    }
                } else {
                    std::cout << "pdgCode: " << pdgCode
                              << " not found in eventID: " << eventID
                              << std::endl;
                    return result;
                }
            }
            ++count;
        }
        std::ranges::shuffle(particlesForThisPDG, generator);
        matrix.resize(particlesForThisPDG.size(), 11);
        for (size_t j = 0; j < particlesForThisPDG.size(); ++j) {
            matrix(j, 0)  = particlesForThisPDG[j].pdg;
            matrix(j, 1)  = particlesForThisPDG[j].px;
            matrix(j, 2)  = particlesForThisPDG[j].py;
            matrix(j, 3)  = particlesForThisPDG[j].pz;
            matrix(j, 4)  = particlesForThisPDG[j].mass;
            matrix(j, 5)  = particlesForThisPDG[j].p0;
            matrix(j, 6)  = particlesForThisPDG[j].x;
            matrix(j, 7)  = particlesForThisPDG[j].y;
            matrix(j, 8)  = particlesForThisPDG[j].z;
            matrix(j, 9)  = particlesForThisPDG[j].freeze_out_time;
            matrix(j, 10) = particlesForThisPDG[j].probability;
        }
        result.push_back(std::move(matrix));
        particlesForThisPDG.clear();
        matrix.resize(0, 11);
    }
    return result;
}
std::tuple<double, double, double, double, double, double>
Coal::fourBodyJacobi(const Coal::Particle &p1, const Coal::Particle &p2,
                     const Coal::Particle &n1, const Coal::Particle &n2) {
    double dx_p1p2  = (p1.x - p2.x) / sqrt(2.0);
    double dy_p1p2  = (p1.y - p2.y) / sqrt(2.0);
    double dz_p1p2  = (p1.z - p2.z) / sqrt(2.0);
    double dpx_p1p2 = (p2.mass * p1.px - p1.mass * p2.px) /
                      (p1.mass + p2.mass) * sqrt(2.0);
    double dpy_p1p2 = (p2.mass * p1.py - p1.mass * p2.py) /
                      (p1.mass + p2.mass) * sqrt(2.0);
    double dpz_p1p2 = (p2.mass * p1.pz - p1.mass * p2.pz) /
                      (p1.mass + p2.mass) * sqrt(2.0);
    double diff_dr_p1p2 =
            sqrt(dx_p1p2 * dx_p1p2 + dy_p1p2 * dy_p1p2 + dz_p1p2 * dz_p1p2);
    double diff_dp_p1p2 = sqrt(dpx_p1p2 * dpx_p1p2 + dpy_p1p2 * dpy_p1p2 +
                               dpz_p1p2 * dpz_p1p2);

    double dx_p1p2n1 =
            (p1.mass * p1.x + p2.mass * p2.x - (p1.mass + p2.mass) * n1.x) /
            (p1.mass + p2.mass) * sqrt(2.0 / 3.0);
    double dy_p1p2n1 =
            (p1.mass * p1.y + p2.mass * p2.y - (p1.mass + p2.mass) * n1.y) /
            (p1.mass + p2.mass) * sqrt(2.0 / 3.0);
    double dz_p1p2n1 =
            (p1.mass * p1.z + p2.mass * p2.z - (p1.mass + p2.mass) * n1.z) /
            (p1.mass + p2.mass) * sqrt(2.0 / 3.0);
    double dpx_p1p2n1 =
            (n1.mass * p1.px + n1.mass * p2.px - (p1.mass + p2.mass) * n1.px) /
            (p1.mass + p2.mass + n1.mass) * sqrt(6) / 2;
    double dpy_p1p2n1 =
            (n1.mass * p1.py + n1.mass * p2.py - (p1.mass + p2.mass) * n1.py) /
            (p1.mass + p2.mass + n1.mass) * sqrt(6) / 2;
    double dpz_p1p2n1 =
            (n1.mass * p1.pz + n1.mass * p2.pz - (p1.mass + p2.mass) * n1.pz) /
            (p1.mass + p2.mass + n1.mass) * sqrt(6) / 2;
    double diff_dr_p1p2n1 = sqrt(dx_p1p2n1 * dx_p1p2n1 + dy_p1p2n1 * dy_p1p2n1 +
                                 dz_p1p2n1 * dz_p1p2n1);
    double diff_dp_p1p2n1 =
            sqrt(dpx_p1p2n1 * dpx_p1p2n1 + dpy_p1p2n1 * dpy_p1p2n1 +
                 dpz_p1p2n1 * dpz_p1p2n1);

    double dx_p1p2n1n2 = (p1.mass * p1.x + p2.mass * p2.x + n1.mass * n1.x -
                          (p1.mass + p2.mass + n1.mass) * n2.x) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dy_p1p2n1n2 = (p1.mass * p1.y + p2.mass * p2.y + n1.mass * n1.y -
                          (p1.mass + p2.mass + n1.mass) * n2.y) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dz_p1p2n1n2 = (p1.mass * p1.z + p2.mass * p2.z + n1.mass * n1.z -
                          (p1.mass + p2.mass + n1.mass) * n2.z) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dpx_p1p2n1n2 = (n2.mass * (p1.px + p2.px + n1.px) -
                           (p1.mass + p2.mass + n1.mass) * n2.px) /
                          (p1.mass + p2.mass + n1.mass + n2.mass) *
                          sqrt(4.0 / 3.0);
    double dpy_p1p2n1n2 = (n2.mass * (p1.py + p2.py + n1.py) -
                           (p1.mass + p2.mass + n1.mass) * n2.py) /
                          (p1.mass + p2.mass + n1.mass + n2.mass) *
                          sqrt(4.0 / 3.0);
    double dpz_p1p2n1n2 = (n2.mass * (p1.pz + p2.pz + n1.pz) -
                           (p1.mass + p2.mass + n1.mass) * n2.pz) /
                          (p1.mass + p2.mass + n1.mass + n2.mass) *
                          sqrt(4.0 / 3.0);
    double diff_dr_p1p2n1n2 =
            sqrt(dx_p1p2n1n2 * dx_p1p2n1n2 + dy_p1p2n1n2 * dy_p1p2n1n2 +
                 dz_p1p2n1n2 * dz_p1p2n1n2);
    double diff_dp_p1p2n1n2 =
            sqrt(dpx_p1p2n1n2 * dpx_p1p2n1n2 + dpy_p1p2n1n2 * dpy_p1p2n1n2 +
                 dpz_p1p2n1n2 * dpz_p1p2n1n2);
    return std::make_tuple(diff_dr_p1p2, diff_dr_p1p2n1, diff_dr_p1p2n1n2,
                           diff_dp_p1p2, diff_dp_p1p2n1, diff_dp_p1p2n1n2);
}

