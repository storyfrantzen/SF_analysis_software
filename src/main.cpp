
#include <iostream>
#include <vector>
#include <string>
#include <cstring>    
#include <algorithm>
#include <filesystem>

#include "includes.h"

namespace fs = std::filesystem;

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
    if (configPath.extension() != ".conf") {
        std::cerr << "ERROR: First argument must be a valid .conf file." << std::endl;
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

    // Must find configuration file containing processing tags:
    std::string CONFIG_FILE = args.configFile;

    try {
        ConfigLoader config(CONFIG_FILE);
    } 
    catch (const std::exception& e) {
        std::cerr << "Failed to load config file: " << e.what() << std::endl;
        return 1;
    }

    // Initialize ConfigLoader:
    ConfigLoader config(CONFIG_FILE);

    // Add numFiles .hipo files to HipoChain:
    clas12root::HipoChain chain;
    for (int i = 0; i < args.numFiles; ++i) { chain.Add(args.files[i]); } 

    // Prepare physics analysis:
    auto config_c12=chain.GetC12Reader();
    auto& c12=chain.C12ref();

    // If input file is NOT MC, populate C12Reader with QA tags found by ConfigLoader:
    if (!config.contains("isMC") || (config.contains("isMC") && !config.getBool("isMC"))) {

        // Populate C12Reader with QA tags found by ConfigLoader:
        if (config_c12->qadb() != nullptr) {
            for (const auto& tag : config.getQARequirements()) {
                config_c12->db()->qadb_addQARequirement(tag);
            }
            //config_c12->db()->qadb_requireOkForAsymmetry(true);
            config_c12->applyQA();
        }
    }

    // Instantiate ProcessManager for workflow:
    ProcessManager PM;

    // Get parameters to filter events:
    const int TORUS = config.getTorus();
    PM.setTorus(TORUS);
    std::string CHANNEL = config.getChannel();
    PM.setChannel(CHANNEL);
    const double EBEAM = config.getEbeam();
    PM.setEbeam(EBEAM);
    std::vector<std::string> TOPOLOGY = config.getTopology();
    PM.setTopology(TOPOLOGY);

    // Initialize output ROOT Tree. This function must be called after CHANNEL_ & TOPOLOGY_ have been set
    std::string OUTPUT_ROOT_FILE = config.getOutputFile();
    PM.rootTree(OUTPUT_ROOT_FILE);

    // Instantiate FiducialCuts, then pass to ProcessManager:
    FiducialCuts FC;

    for (const auto& tag : config.getFiducialCuts()) { FC.addCut(tag); }
    PM.setFiducialCuts(FC);

    int MAX_EVENTS = 10000000;
    int counter = 0;
    while (chain.Next() && counter < MAX_EVENTS) { 
        PM.processEvent(*c12);
        counter++;
    }
    PM.finalize();

    if (config.contains("writeSummary") && config.getBool("writeSummary")) config.writeSummary();

    std::cout << "Processing complete." << std::endl;
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;

    return 0;
}