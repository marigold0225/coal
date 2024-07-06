//
// Created by mafu on 1/8/2024.
//
#include "../headers/method.h"
#include "../headers/random.h"
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <ranges>

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
    if (label == "0-100") {
        return filename + "/" + fileType + ".dat";
    }
    return filename + "/" + label + "/" + fileType + "_" + label + ".dat";
}

std::vector<double> Coal::linspace(const double start, const double end, const int num) {
    std::vector<double> Linspaced;
    const double delta = (end - start) / (num - 1);

    Linspaced.reserve(num - 1);
    for (int i = 0; i < num - 1; ++i) {
        Linspaced.push_back(start + delta * i);
    }
    Linspaced.push_back(end);

    return Linspaced;
}

void Coal::lorenzBoostMatrix(Eigen::Matrix<double, Eigen::Dynamic, 11> &particles,
                             const double beta_x, const double beta_y,
                             const double beta_z) {
    double beta2                = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;
    if (beta2 == 0 || beta2 < 1e-6) {
        return;
    }
    if (beta2 > 0.999999999) {
        beta2 = 0.999999999;
    }
    const double gamma          = 1.0 / sqrt(1.0 - beta2);
    const double _gamma_1_beta2 = (gamma - 1) / beta2;
    Eigen::Matrix4d lambda;
    lambda << gamma, -gamma * beta_x, -gamma * beta_y, -gamma * beta_z, -gamma * beta_x,
            1 + _gamma_1_beta2 * beta_x * beta_x, _gamma_1_beta2 * beta_x * beta_y,
            _gamma_1_beta2 * beta_x * beta_z, -gamma * beta_y,
            _gamma_1_beta2 * beta_y * beta_x, 1 + _gamma_1_beta2 * beta_y * beta_y,
            _gamma_1_beta2 * beta_y * beta_z, -gamma * beta_z,
            _gamma_1_beta2 * beta_z * beta_x, _gamma_1_beta2 * beta_z * beta_y,
            1 + _gamma_1_beta2 * beta_z * beta_z;
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
void Coal::boostToComMatrix(Eigen::Matrix<double, Eigen::Dynamic, 11> &particles) {
    const double beta_x = particles.col(1).sum() / particles.col(5).sum();
    const double beta_y = particles.col(2).sum() / particles.col(5).sum();
    const double beta_z = particles.col(3).sum() / particles.col(5).sum();

    lorenzBoostMatrix(particles, beta_x, beta_y, beta_z);

    const double t_max = particles.col(9).maxCoeff();

    for (int i = 0; i < particles.rows(); ++i) {
        particles(i, 6) += (t_max - particles(i, 9)) * particles(i, 1) / particles(i, 5);
        particles(i, 7) += (t_max - particles(i, 9)) * particles(i, 2) / particles(i, 5);
        particles(i, 8) += (t_max - particles(i, 9)) * particles(i, 3) / particles(i, 5);
    }
}

void Coal::calculateLorentz(
        Eigen::Ref<Eigen::VectorXd> betaX, Eigen::Ref<Eigen::VectorXd> betaY,
        Eigen::Ref<Eigen::VectorXd> betaZ,
        std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
                &lorentzMatrixs) {

    assert(lorentzMatrixs.size() == betaX.size());
    for (int i = 0; i < betaX.size(); ++i) {
        const double beta_x = betaX(i);
        const double beta_y = betaY(i);
        const double beta_z = betaZ(i);
        double beta2        = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;

        if (beta2 == 0 || beta2 < 1e-6) {
            lorentzMatrixs[i] = Eigen::Matrix4d::Identity();
            continue;
        }
        if (beta2 > 0.999999999) {
            beta2 = 0.999999999;
        }
        double gamma                = 1.0 / sqrt(1.0 - beta2);
        const double _gamma_1_beta2 = (gamma - 1) / beta2;
        Eigen::Matrix4d lambda;
        lambda << gamma, -gamma * beta_x, -gamma * beta_y, -gamma * beta_z,
                -gamma * beta_x, 1 + _gamma_1_beta2 * beta_x * beta_x,
                _gamma_1_beta2 * beta_x * beta_y, _gamma_1_beta2 * beta_x * beta_z,
                -gamma * beta_y, _gamma_1_beta2 * beta_y * beta_x,
                1 + _gamma_1_beta2 * beta_y * beta_y, _gamma_1_beta2 * beta_y * beta_z,
                -gamma * beta_z, _gamma_1_beta2 * beta_z * beta_x,
                _gamma_1_beta2 * beta_z * beta_y, 1 + _gamma_1_beta2 * beta_z * beta_z;
        lorentzMatrixs[i] = lambda;
    }
}

void Coal::applyLorentzBoost(
        Eigen::Ref<Eigen::MatrixXd> combinedX, Eigen::Ref<Eigen::MatrixXd> combinedY,
        Eigen::Ref<Eigen::MatrixXd> combinedZ, Eigen::Ref<Eigen::MatrixXd> combinedT,
        Eigen::Ref<Eigen::MatrixXd> combinedPX, Eigen::Ref<Eigen::MatrixXd> combinedPY,
        Eigen::Ref<Eigen::MatrixXd> combinedPZ, Eigen::Ref<Eigen::MatrixXd> combinedP0,
        const std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
                &lorentz) {

    Eigen::Vector4d particleMomentum, particlePosition, boostedMomentum, boostedPosition;

    for (int j = 0; j < combinedX.cols(); ++j) {
        for (int i = 0; i < combinedX.rows(); ++i) {


            particleMomentum << combinedP0(i, j), combinedPX(i, j), combinedPY(i, j),
                    combinedPZ(i, j);
            particlePosition << combinedT(i, j), combinedX(i, j), combinedY(i, j),
                    combinedZ(i, j);

            boostedMomentum.noalias() = lorentz[i] * particleMomentum;
            boostedPosition.noalias() = lorentz[i] * particlePosition;

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
void Coal::applyLorentzBoost_col(
        Eigen::Ref<Eigen::MatrixXd> combinedX, Eigen::Ref<Eigen::MatrixXd> combinedY,
        Eigen::Ref<Eigen::MatrixXd> combinedZ, Eigen::Ref<Eigen::MatrixXd> combinedT,
        Eigen::Ref<Eigen::MatrixXd> combinedPX, Eigen::Ref<Eigen::MatrixXd> combinedPY,
        Eigen::Ref<Eigen::MatrixXd> combinedPZ, Eigen::Ref<Eigen::MatrixXd> combinedP0,
        const std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
                &lorentz) {
    Eigen::Vector4d particleMomentum, particlePosition, boostedMomentum, boostedPosition;
    assert(combinedX.rows() == 4 * lorentz.size());
    for (int j = 0; j < combinedX.cols(); ++j) {
        for (int i = 0; i < combinedX.rows(); ++i) {

            particleMomentum << combinedP0(i, j), combinedPX(i, j), combinedPY(i, j),
                    combinedPZ(i, j);
            particlePosition << combinedT(i, j), combinedX(i, j), combinedY(i, j),
                    combinedZ(i, j);

            boostedMomentum.noalias() = lorentz[i] * particleMomentum;
            boostedPosition.noalias() = lorentz[i] * particlePosition;

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

std::tuple<std::vector<double>, std::vector<double>>
Coal::JacobiCoordinatesMatrix_test(const Eigen::MatrixXd &particles,
                                   const ClusterParams &params) {
    const auto N = particles.rows();
    Eigen::MatrixXd M, M_inv_t;
    if (!params.Fixed) {
        std::vector<double> massArray(N);
        for (int i = 0; i < N; ++i) {
            massArray[i] = particles(i, 4);
        }
        M       = MassMatrix(massArray);
        M_inv_t = M.inverse().transpose();
    } else {
        if (N - 2 < params.M.size()) {
            M       = params.M[N - 2];
            M_inv_t = params.M_inv_t[N - 2];
        } else {
            throw std::runtime_error("Particle count exceeds precomputed matrix size.");
        }
    }
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

    Eigen::VectorXd x_jacobi  = transformCoordinates(x, M);
    Eigen::VectorXd y_jacobi  = transformCoordinates(y, M);
    Eigen::VectorXd z_jacobi  = transformCoordinates(z, M);
    Eigen::VectorXd px_jacobi = transformCoordinates(px, M_inv_t);
    Eigen::VectorXd py_jacobi = transformCoordinates(py, M_inv_t);
    Eigen::VectorXd pz_jacobi = transformCoordinates(pz, M_inv_t);

    std::vector<double> d_r(N - 1), d_p(N - 1);
    for (auto i = 0; i < N - 1; ++i) {
        d_r[i] = sqrt(x_jacobi(i) * x_jacobi(i) + y_jacobi(i) * y_jacobi(i) +
                      z_jacobi(i) * z_jacobi(i));
        d_p[i] = sqrt(px_jacobi(i) * px_jacobi(i) + py_jacobi(i) * py_jacobi(i) +
                      pz_jacobi(i) * pz_jacobi(i));
    }
    return std::make_tuple(d_r, d_p);
}

std::tuple<std::vector<double>, std::vector<double>>
Coal::JacobiCoordinatesMatrix(const Eigen::MatrixXd &particles,
                              const ClusterParams &params) {
    const auto N = particles.rows();

    auto transformCoordinates = [&](const Eigen::VectorXd &vec,
                                    const Eigen::MatrixXd &mat) {
        return (mat * vec).tail(vec.size() - 1);
    };

    Eigen::Matrix<double, Eigen::Dynamic, 3> r_jacobi(N - 1, 3), p_jacobi(N - 1, 3);
    r_jacobi.col(0) = transformCoordinates(particles.col(6), params.M[N - 2]);
    r_jacobi.col(1) = transformCoordinates(particles.col(7), params.M[N - 2]);
    r_jacobi.col(2) = transformCoordinates(particles.col(8), params.M[N - 2]);
    p_jacobi.col(0) = transformCoordinates(particles.col(1), params.M_inv_t[N - 2]);
    p_jacobi.col(1) = transformCoordinates(particles.col(2), params.M_inv_t[N - 2]);
    p_jacobi.col(2) = transformCoordinates(particles.col(3), params.M_inv_t[N - 2]);

    std::vector<double> d_r(N - 1), d_p(N - 1);

    for (int i = 0; i < N - 1; ++i) {
        d_r[i] = r_jacobi.row(i).norm();
        d_p[i] = p_jacobi.row(i).norm();
    }

    return {d_r, d_p};
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd>
Coal::JacobiCoordinatesMatrix_Vec(const Eigen::MatrixXd &particles,
                                  const ClusterParams &params) {
    const auto N = particles.rows();

    Eigen::Matrix<double, Eigen::Dynamic, 3> r_jacobi(N - 1, 3), p_jacobi(N - 1, 3);
    r_jacobi.col(0) = (params.M[N - 2] * particles.col(6)).tail(N - 1);
    r_jacobi.col(1) = (params.M[N - 2] * particles.col(7)).tail(N - 1);
    r_jacobi.col(2) = (params.M[N - 2] * particles.col(8)).tail(N - 1);
    p_jacobi.col(0) = (params.M_inv_t[N - 2] * particles.col(1)).tail(N - 1);
    p_jacobi.col(1) = (params.M_inv_t[N - 2] * particles.col(2)).tail(N - 1);
    p_jacobi.col(2) = (params.M_inv_t[N - 2] * particles.col(3)).tail(N - 1);

    Eigen::VectorXd d_r = r_jacobi.rowwise().norm();
    Eigen::VectorXd d_p = p_jacobi.rowwise().norm();

    return {d_r, d_p};
}

Coal::EventsMap Coal::selectEvents(const EventsMap &eventMap, const ClusterParams &params,
                                   ResParamsMap &resolution, const int currentIndex,
                                   std::vector<int> &eventIDList) {
    EventsMap sampleEvents;
    if (params.mixEvents == 1) {
        const int selectIndex       = currentIndex % static_cast<int>(eventMap.size());
        const int selectEventId     = eventIDList[selectIndex];
        sampleEvents[selectEventId] = eventMap.at(selectEventId);
        resolution.selectEventID[currentIndex] = {selectEventId, 0};

    } else {

        int selectSize  = params.mixEvents * params.NBody;
        auto &generator = RandomNumber::getInstance().getGenerator();
        std::ranges::shuffle(eventIDList, generator);
        selectSize = std::min(selectSize, static_cast<int>(eventIDList.size()));
        for (int i = 0; i < selectSize; ++i) {
            auto eventID          = eventIDList[i];
            sampleEvents[eventID] = eventMap.at(eventID);
            const double psi_1    = resolution.eventPlaneMap.at(eventID);
            rotateEventPlane(sampleEvents[eventID], -psi_1);
        }
    }
    return sampleEvents;
}

Coal::EventsMap Coal::selectEvents_v2(const EventsMap &eventMap,
                                      const ClusterParams &params,
                                      std::vector<int> &eventIDList) {
    EventsMap sampleEvents;

    int selectSize  = params.mixEvents * params.NBody;
    auto &generator = RandomNumber::getInstance().getGenerator();
    std::ranges::shuffle(eventIDList, generator);
    selectSize = std::min(selectSize, static_cast<int>(eventIDList.size()));
    for (int i = 0; i < selectSize; ++i) {
        auto eventID          = eventIDList[i];
        sampleEvents[eventID] = eventMap.at(eventID);
    }
    return sampleEvents;
}

std::vector<Eigen::MatrixXd> Coal::selectParticles(const EventsMap &eventMap,
                                                   const ClusterParams &params) {
    const auto PDGs = params.PDGArray;
    std::vector<Eigen::MatrixXd> result;
    auto &generator = RandomNumber::getInstance().getGenerator();

    if (eventMap.size() == 1) {
        const auto &event = eventMap.begin()->second;
        for (auto i = 0; i < PDGs.size(); ++i) {
            const auto pdgCode = PDGs[i];
            ParticleArray particlesForThisPDG;
            if (auto particleIt = event.find(pdgCode); particleIt != event.end()) {
                for (const auto &particle: particleIt->second) {
                    if (const double rapidity = particle.getRapidity();
                        rapidity > params.originRapidity[i].first &&
                        rapidity < params.originRapidity[i].second) {
                        particlesForThisPDG.push_back(particle);
                    }
                }
            } else {
                std::cout << "pdgCode: " << pdgCode
                          << " not found in eventID: " << eventMap.begin()->first
                          << std::endl;
                continue;
            }
            std::ranges::shuffle(particlesForThisPDG, generator);
            Eigen::MatrixXd matrix(static_cast<long>(particlesForThisPDG.size()), 11);
            for (long j = 0; j < particlesForThisPDG.size(); ++j) {
                const auto &p = particlesForThisPDG[j];
                matrix.row(j) = p.convertToMatrix();
            }
            result.push_back(std::move(matrix));
        }
    } else {
        std::vector<int> eventIDList;
        std::ranges::transform(eventMap, std::back_inserter(eventIDList),
                               [](const auto &pair) { return pair.first; });
        std::ranges::shuffle(eventIDList, generator);
        ParticleArray particlesForThisPDG{};
        Eigen::MatrixXd matrix(0, 11);
        for (size_t i = 0; i < PDGs.size(); ++i) {
            int pdgCode = PDGs[i];
            int count   = 0;

            while (count < params.mixEvents && !eventIDList.empty()) {
                int eventID = eventIDList.back();
                eventIDList.pop_back();

                if (auto eventIt = eventMap.find(eventID); eventIt != eventMap.end()) {
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
                                  << " not found in eventID: " << eventID << std::endl;
                        continue;
                    }
                }
                ++count;
            }
            std::ranges::shuffle(particlesForThisPDG, generator);
            matrix.resize(static_cast<long>(particlesForThisPDG.size()), 11);
            for (long j = 0; j < particlesForThisPDG.size(); ++j) {
                const auto &p = particlesForThisPDG[j];
                matrix.row(j) = p.convertToMatrix();
            }
            result.push_back(std::move(matrix));
            particlesForThisPDG.clear();
            matrix.resize(0, 11);
        }
    }
    return result;
}
