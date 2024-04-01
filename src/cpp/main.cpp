//
// Created by mafu on 1/2/2024.
//
#include "../headers/run.h"
#include <filesystem>
#include <fstream>
#include <iostream>

int main(const int argc, char* argv[]) {
    std::cout << "currentPath is: " << std::filesystem::current_path() << std::endl;

    std::string inputfile = "input/input.yaml";
    for (int i = 1; i < argc - 1; ++i) {
        if (std::string(argv[i]) == "-i") {
            inputfile = argv[i + 1];
        }
    }

    std::cout << "Using input file: " << inputfile << std::endl;

    auto config = Coal::ConfigParser(inputfile);
    const auto allEvents = Coal::PreData::readFile(
        config.inputFileName, config.general["Mode"].as<std::string>());
    handleReactions(allEvents, config);

    std::cout << "All calculating finished!" << std::endl;

    return 0;
}

// int main() {
//
//     std::cout << "currentPath is:" << std::filesystem::current_path()
//               << std::endl;
//
//     const std::string inputfile = "input/input.yaml";
//     auto config          = Coal::ConfigParser(inputfile);
//     const auto allEvents = Coal::PreData::readFile(
//             config.inputFileName, config.general["Mode"].as<std::string>());
//     handleReactions(allEvents, config);
//     std::cout << "All calculating finshed!" << std::endl;
//     return 0;
// }
