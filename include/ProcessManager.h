#pragma once

#include "nlohmann/json.hpp"
#include "FiducialCuts.h"  // Note: ProcessManager has its own FC object during filtering
#include "BranchVars.h" // IMPORTANT: contains structs that manage ALL quantities to be stored
#include "PhysicalConstants.h" // contains useful constants
#include "clas12reader.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TObject.h"
#include <string>
#include <cstdint>  // for int8_t
#include <iomanip>  // for std::put_time
#include <sstream>  // for std::ostringstream
#include <ctime>    // for std::time_t, std::localtime

using namespace clas12;

class ProcessManager {
public:

    ProcessManager(const nlohmann::json& config);

    int eventsProcessed() {return eventsProcessed_;}
    int numFills() {return numFills_;}

    std::string currentTimestamp() const; 
    static int getDetector(int status);
    bool passesVertexCut(const float vz, const int zmin=-8, const int zmax=2);

    // saves ROOT tree to output file: 
    void finalize();

    // BIGGER FUNCTIONS: //
    void processEvent(clas12::clas12reader& c12);
    void processEPPI0(clas12::clas12reader& c12);

private:

    int eventsProcessed_ = 0;
    int numFills_ = 0;

    double ebeam_;
    std::string channel_;
    int torus_;
    std::string outPrefix_;
    bool requireTopology_ = false;

    FiducialCuts* FC_ = nullptr;
    TTree* tree_      = nullptr;

    EventVars ev_;
    RecVars   rec_;
    RecVars   e_;
    RecVars   p_;
    RecVars   g_;
    GenVars   gen_;
    DISVars   dis_;
    DISVars   gen_dis_;
    EPPI0Vars eppi0_;

    std::unordered_set<std::string> enabledEvBranches_  = {};
    std::unordered_set<std::string> enabledRecBranches_ = {};
    std::unordered_set<std::string> enabledEleBranches_ = {};
    std::unordered_set<std::string> enabledProBranches_ = {};
    std::unordered_set<std::string> enabledPhoBranches_ = {};

};

