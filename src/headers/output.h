//
// Created by mafu on 4/2/2024.
//
#pragma once
#include "flow.h"
#include "shortcut.h"
#include <ranges>

#include <fstream>
#include <iostream>

namespace Coal {

    class outputInterface final {
    public:
        explicit outputInterface(const std::string &filename)
            : is_pt_output(filename == "pt"), is_flow_output(filename == "flow"),
              is_cluster_output(filename == "cluster") {}

        ~outputInterface() = default;

        void output(const Eigen::MatrixXd &targetParticles,
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
                   << result.v1overpt[rap] / result.v1overptCount[rap]// / resolution.Ref_v1
                   << ", v2:"
                   << result.v2overpt[rap] / result.v2overptCount[rap]// / resolution.Ref_v2
                   << "\n";
            for (int i = 0; i < result.ptBins.second; ++i) {
                const double pt = result.ptBins.first / 2 + i * result.ptBins.first;
                result.ptArray[rap][i] /= (2 * M_PI * pt * result.ptBins.first *
                                           std::abs(rap.second - rap.first));

                result.v1Array[rap][i] = result.v1Count[rap][i] > 0
                                                 ? result.v1Array[rap][i] /
                                                           result.v1Count[rap][i]
                // /resolution.Ref_v1
                                                 : 0.0;
                result.v2Array[rap][i] = result.v2Count[rap][i] > 0
                                                 ? result.v2Array[rap][i] /
                                                           result.v2Count[rap][i]
                // /resolution.Ref_v2
                                                 : 0.0;

                output << pt << " " << result.ptArray[rap][i] << " "
                       << result.v1Array[rap][i] << " " << result.v2Array[rap][i] << "\n";
            }
        }
    };
}// namespace Coal
