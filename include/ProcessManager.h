#pragma once

#include "nlohmann/json.hpp"
#include "FiducialCuts.h"  // Note: ProcessManager has its own FC object during filtering
#include "KinematicCorrections.h" // ProcessManager has a member instance of KC during filtering
#include "SFCuts.h" // ProcessManager has a member instance of SF during filtering
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
    bool channelCheck(float Q2, float W, float y);
    bool passesVertexCut(const int pid, const float vz, const float e_vz);

    // saves ROOT tree to output file and creates summary tree which stores Events tree and Summary tree: 
    void finalize(double totalCharge);

    // BIGGER FUNCTIONS: //
    void processEvent(clas12::clas12reader& c12);
    void processEPPI0REC(clas12::clas12reader& c12);
    void processEPPI0GEMC(clas12::clas12reader& c12);
    void processEPPI0MATCH(clas12::clas12reader& c12);

private:

    int eventsProcessed_ = 0;
    int numFills_ = 0;

    double ebeam_;
    std::string channel_;
    int torus_;
    std::string tag_;
    bool requireTopology_ = false;

    std::unique_ptr<FiducialCuts> FC_;
    std::unique_ptr<KinematicCorrections> KC_;
    std::unique_ptr<SFCuts> SF_;
    
    TFile* outFile_   = nullptr;
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
    EPPI0Vars gen_eppi0_;

    std::unordered_set<std::string> enabledEvBranches_  = {};
    std::unordered_set<std::string> enabledRecBranches_ = {};
    std::unordered_set<std::string> enabledEleBranches_ = {};
    std::unordered_set<std::string> enabledProBranches_ = {};
    std::unordered_set<std::string> enabledPhoBranches_ = {};

};

