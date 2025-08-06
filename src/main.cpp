
#include <iostream>
#include <vector>
#include <string>
#include <cstring>    
#include <algorithm>
#include <filesystem>

#include "includes.h"

#include "nlohmann/json.hpp"
#include <fstream>

namespace fs = std::filesystem;
using json = nlohmann::json;
using namespace std;


struct ArgResults {
    std::vector<std::string> files;
    std::string configFile;
    int numFiles;
};
ArgResults readArgs(int argc, char** argv) {
    ArgResults result;

    if (argc < 3) {
        std::cerr << "Usage (from project build directory): ./processHipoFiles <config.conf> <hipo_dir1> [<hipo_dir2> ...] [<max_files:int>]" << std::endl;
        return result; // returns empty files and numFiles = 0
    }

    auto is_number = [](const char* s) -> bool {
        return s && *s && std::all_of(s, s + std::strlen(s), ::isdigit);
    };

    result.numFiles = -1;
    int numDirs = argc - 2;

    // Determine if the last argument is a number: if not, numDirs = arc - 2; if so, numDirs = argc - 3.
    if (is_number(argv[argc - 1])) {
        result.numFiles = std::stoi(argv[argc - 1]);
        numDirs--;
    }
    if (numDirs < 1) {
        std::cerr << "ERROR: You must provide at least one hipo directory." << std::endl;
        return result;
    }
    fs::path configPath(argv[1]);
    if (configPath.extension() != ".json") {
        std::cerr << "ERROR: First argument must be a valid .json file." << std::endl;
        return result;
    }
    result.configFile = argv[1];

    for (int i = 2; i < numDirs + 2; ++i) {
        fs::path dir(argv[i]);
        if (!fs::exists(dir) || !fs::is_directory(dir)) {
            std::cerr << "WARNING: " << dir << " is not a valid directory, skipping." << std::endl;
            continue;
        }
        for (const auto& entry : fs::recursive_directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                result.files.push_back(entry.path().string());
            }
        }
    }

    if (result.files.empty()) {
        std::cerr << "ERROR: No .hipo files found in the specified directories." << std::endl;
        return result;
    }

    if (result.numFiles < 0 || result.numFiles > static_cast<int>(result.files.size())) {
        result.numFiles = static_cast<int>(result.files.size());
    }

    std::cout << "Using " << result.numFiles << " files out of " << result.files.size() << " available." << std::endl;

    return result;
}

int main(int argc, char** argv) {

    auto start = std::chrono::high_resolution_clock::now();
    auto args = readArgs(argc, argv);

    // Must find configuration .json file containing processing tags:
    std::ifstream configStream(args.configFile);
    if (!configStream.is_open()) {
        std::cerr << "Failed to open config file: " << args.configFile << std::endl;
        return 1;
    }
    
    json config;
    configStream >> config;

    // Add numFiles .hipo files to HipoChain:
    clas12root::HipoChain chain;
    for (int i = 0; i < args.numFiles; ++i) {
        chain.Add(args.files[i]);
    } 

    // Prepare physics analysis:
    auto config_c12=chain.GetC12Reader();
    auto& c12=chain.C12ref();

    // QA:
    if (!config.value("isMC", false)) {
        if (config_c12->qadb() != nullptr) {
            for (const auto& tag : config["qa"])
                config_c12->db()->qadb_addQARequirement(tag);
            config_c12->applyQA();
        }
    }

    // Instantiate ProcessManager for workflow:
    ProcessManager PM(config);

    // Initialize ROOT Tree. This function must be called after CHANNEL_ & TOPOLOGY_ have been set 
    PM.rootTree();

    int MAX_EVENTS = 10000000;
    while (chain.Next() && PM.eventsProcessed() < MAX_EVENTS) { 
        PM.processEvent(*c12);
        if (PM.eventsProcessed() % 500000 == 0) std::cout << PM.eventsProcessed() << " events processed. \n";
    }

    std::string outFile = config.value("outFile", PM.makeFilename());
    PM.finalize(outFile);

    std::cout << "Processing complete. " << PM.eventsProcessed() << " events were processed, and " <<
                    PM.numFills() << " tree fills were made." << std::endl;
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;

    return 0;
}