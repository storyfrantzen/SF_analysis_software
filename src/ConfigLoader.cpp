#include <iostream>
#include "ConfigLoader.h"

std::vector<std::string> ConfigLoader::getQARequirements() const {
    if (!contains("qadbRequirements")) {
        // Default to no requirements
        return {};
    }

    std::vector<std::string> requirements;
    std::istringstream iss(get("qadbRequirements"));
    std::string token;
    while (std::getline(iss, token, ',')) {
        //std::cout << "Token before trimming : " << token << std::endl;
        trim(token);
        //std::cout << "Token after trimming : " << "[" << token << "]" << std::endl;  // Debug print
        if (!token.empty()) requirements.push_back(token);
    }
    return requirements;
}

std::string ConfigLoader::getOutputFile() const {
    if (contains("output")) return get("output");
    return "/work/clas12/storyf/SF_analysis_software/output/dummy.root";  // DEFAULT output filepath
}

std::string ConfigLoader::getChannel() const {
    if (contains("channel")) return get("channel");
    return "eppi0";  // DEFAULT channel
}

int ConfigLoader::getTorus() const {
    if (contains("torus")) return getInt("torus");
    return +1;  // DEFAULT torus = outbending
}

double ConfigLoader::getEbeam() const {
    if (contains("Ebeam")) return getDouble("Ebeam");
    return 6.535; // DEFAULT beam energy
}

std::vector<std::string> ConfigLoader::getFiducialCuts() const {
    std::vector<std::string> cuts;
    if (!contains("fiducialCuts")) return cuts;

    std::istringstream iss(get("fiducialCuts"));
    std::string token;
    while (std::getline(iss, token, ',')) {
        trim(token); 
        if (!token.empty())
            cuts.push_back(token);
    }
    return cuts;
}

std::vector<std::string> ConfigLoader::getTopology() const {

    std::vector<std::string> topo;
    if (!contains("topology")) return topo;

    std::istringstream iss(get("topology"));
    std::string token;
    while (std::getline(iss, token, ',')) {
        trim(token);
        if (!token.empty())
            topo.push_back(token);
    }
    return topo;
}
