//
// Created by mafu on 1/2/2024.
//
#include "../headers/run.h"
#include <filesystem>
#include <fstream>
#include <iostream>

int main() {

    std::cout << "currentPath is:" << std::filesystem::current_path()
              << std::endl;

    const std::string inputfile = "input/input.yaml";

    //
    // YAML::Node config = YAML::LoadFile(inputfile);
    // if (auto General = config["General"];
    //     General["Centrality"]["Enable"].as<bool>()) {
    //     Coal::handleCentralityCalculation(inputfile);
    // } else {
    //     Coal::handleNoCentralityCalculation(inputfile);
    // }

    auto config          = Coal::readCondig(inputfile);
    const auto allEvents = Coal::PreData::readFile(
            config.inputFileName, config.general["Mode"].as<std::string>());
    handleReactions(allEvents, config);
    std::cout << "All calculating finshed!" << std::endl;
    return 0;
}
