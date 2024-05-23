//
// Created by mafu on 4/2/2024.
//
#pragma once
#include "flow.h"
#include "shortcut.h"
#include <ranges>

#include <fstream>
#include <iostream>
#include <numeric>

namespace Coal {

    class outputInterface final {
    public:
        explicit outputInterface(const std::string &filename)
            : is_pt_output(filename == "pt"), is_flow_output(filename == "flow"),
              is_particle_output(filename == "particle"),
              is_cluster_output(filename == "cluster") {}

        ~outputInterface() = default;

        template<typename MatrixType>
        void output(const Eigen::MatrixBase<MatrixType> &targetParticles,
                    const std::string &outputFileName,
                    const ClusterParams &params) const {
            if (!is_cluster_output) {
                return;
            }

            std::ofstream output(outputFileName);
            if (!output.is_open()) {
                std::cerr << "Error: unable to open file " << outputFileName << std::endl;
                return;
            }
            const double totalEvents = params.eventFactor;
            output << "#events:" << totalEvents << " "
                   << "size:" << targetParticles.rows() << "\n";
            output << targetParticles << "\n";
            output.close();
        }

        void output(const EventsMap &allevents, const std::string &outputFIleName,
                    const int pdg) const {
            if (!is_particle_output) {
                return;
            }
            std::ofstream output(outputFIleName);
            if (!output.is_open()) {
                std::cerr << "Error: unable to open file " << outputFIleName << std::endl;
                return;
            }
            const size_t totalEvents = allevents.size();
            const int total_particles =
                    std::accumulate(allevents.begin(), allevents.end(), 0,
                                    [](const size_t &sum, const auto &event) {
                                        return sum + event.second.size();
                                    });
            output << "#events:" << totalEvents << " "
                   << "size:" << total_particles << "\n";
            for (const auto &OneEvent: allevents | std::ranges::views::values) {
                for (const auto &[particleType, particles]: OneEvent) {
                    if (particleType != pdg) {
                        continue;
                    }
                    for (const auto &particle: particles) {
                        output << particle << "\n";
                    }
                }
            }
            output.close();
        }

        void output(flowResult &result, const std::string &outputFileName,
                    const ResParamsMap &resolution) const {
            if (!is_pt_output && !is_flow_output) {
                return;
            }

            std::ofstream output(outputFileName);
            if (!output.is_open()) {
                std::cerr << "Error: unable to open file " << outputFileName << std::endl;
                return;
            }

            for (const auto &key: result.yeildArray | std::views::keys) {
                const auto &rap = key;
                outputRapidityData(output, result, rap, resolution);
            }

            output.close();
        }

    private:
        const bool is_pt_output;
        const bool is_flow_output;
        const bool is_particle_output;
        const bool is_cluster_output;

        void outputRapidityData(std::ofstream &output, flowResult &result,
                                const RapidityRange &rap,
                                const ResParamsMap &resolution) const {
            output << "Rapidity range: " << rap.first << "<y<" << rap.second
                   << ", cluster yield:"
                   << result.yeildArray[rap] / (rap.second - rap.first) << "\n";

            if (is_flow_output) {
                outputFlowData(output, result, rap, resolution);
            }

            if (is_pt_output) {
                outputPTData(output, result, rap, resolution);
            }

            output << "\n";
        }

        static void outputFlowData(std::ofstream &output, flowResult &result,
                                   const RapidityRange &rap,
                                   const ResParamsMap &resolution) {
            output << "flow over pt range: " << result.ptRange.first << "<pt<"
                   << result.ptRange.second << ", v1:"
                   << result.v1overpt[rap] / result.v1overptCount[rap] / resolution.Ref_v1
                   << ", v2:"
                   << result.v2overpt[rap] / result.v2overptCount[rap] / resolution.Ref_v2
                   << "\n";
        }

        static void outputPTData(std::ofstream &output, flowResult &result,
                                 const RapidityRange &rap,
                                 const ResParamsMap &resolution) {
            output << "flow over pt range: " << result.ptRange.first << "<pt<"
                   << result.ptRange.second << ", v1:"
                   << result.v1overpt[rap] / result.v1overptCount[rap] / resolution.Ref_v1
                   << ", v2:"
                   << result.v2overpt[rap] / result.v2overptCount[rap] / resolution.Ref_v2
                   << "\n";
            for (int i = 0; i < result.ptBins.second; ++i) {
                const double pt = result.ptBins.first / 2 + i * result.ptBins.first;
                // result.ptArray[rap][i] /= (2 * M_PI * pt * result.ptBins.first *
                //                            std::abs(rap.second - rap.first));
                result.ptArray[rap][i] /= ( result.ptBins.first *
                                           std::abs(rap.second - rap.first));
                result.v1Array[rap][i] = (result.v1Count[rap][i] > 0 && resolution.Ref_v1 != 0.0)
                                                 ? result.v1Array[rap][i] /
                                                           result.v1Count[rap][i] /
                                                           resolution.Ref_v1
                                                 : 0.0;
                result.v2Array[rap][i] = (result.v2Count[rap][i] > 0 && resolution.Ref_v2 != 0.0)
                                                 ? result.v2Array[rap][i] /
                                                           result.v2Count[rap][i] /
                                                           resolution.Ref_v2
                                                 : 0.0;

                output << pt << " " << result.ptArray[rap][i] << " "
                       << result.v1Array[rap][i] << " " << result.v2Array[rap][i] << "\n";
            }
        }
    };
}// namespace Coal
