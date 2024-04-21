//
// Created by mafu on 1/2/2024.
//
#include "../headers/logger.h"
#include "../headers/run.h"
#include <filesystem>
#include <iostream>

int main(const int argc, char *argv[]) {

    std::string inputfile = "input/input.yaml";
    for (int i = 1; i < argc - 1; ++i) {
        if (std::string(argv[i]) == "-i") {
            inputfile = argv[i + 1];
        }
    }

    if (!std::filesystem::exists(inputfile)) {
        std::cerr << "Input file does not exist: " << inputfile << std::endl;
        return 1;
    }

    const auto config          = Coal::ConfigParser(inputfile);

    const auto loggerPath = config.outputPath + "/logfile.log";

    Coal::Logger::getInstance().init(loggerPath);

    const auto logger = Coal::Logger::getInstance().get();

    logger->info("Current path is: {}", std::filesystem::current_path().string());

    const auto allEvents = Coal::FileLoader::readFile(config.inputFileName, config.mode);

    handleReactions(allEvents, config);

    logger->info("All calculating finished!");

    return 0;
}
