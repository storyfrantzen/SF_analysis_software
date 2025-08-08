#pragma once

#include "nlohmann/json.hpp"
#include "FiducialCuts.h"  // Note: ProcessManager has its own FC object during filtering
#include "Kinematics.h" // Note: ProcessManager will use kinematics class for computing derived quantities
#include "Vars.h"
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

    // Getters defined here: //
    int eventsProcessed() {return eventsProcessed_;}
    int numFills() {return numFills_;}
    // Based on particle status, returns int representation of FT (0), FD (1), or CD (2)
    static int getDetector(int status) {
        const int absStatus = std::abs(status);

        if (absStatus >= 1000 && absStatus < 2000) return 0; // FT
        if (absStatus >= 2000 && absStatus < 4000) return 1; // FD
        if (absStatus >= 4000 && absStatus < 5000) return 2; // CD

        return -999; // Unknown detector
    }

    // LITTLE FUNCTIONS: //
    std::string currentTimestamp() const; 
    std::string makeFilename() const;
    bool passesVertexCut(clas12::region_particle* p, const int zmin=-8, const int zmax=2);
    bool passesDiagECALCut(clas12::region_particle* ele);

    // BIGGER FUNCTIONS: //
    void writeEleBranches(int info=1);
    void writeProBranches(int info=1);
    void writePhoBranches(int info=1);
    void writeEPPI0Branches();
    void rootTree();
    void fillRecVars(clas12::region_particle* p, int ele_info=1, int pro_info=1, int pho_info=1);
    void fillEleVars(clas12::region_particle* ele, int info=1);
    void fillProVars(clas12::region_particle* pro, int info=1);
    void fillPhoVars(clas12::region_particle* pho, int info=1);
    void fillEPPI0Vars(Kinematics EPPI0);
    void processEPPI0(clas12::clas12reader& c12);
    void processEvent(clas12::clas12reader& c12);
    void finalize(const std::string& output_file);

private:

    int eventsProcessed_ = 0;
    int numFills_ = 0;

    double ebeam_;
    std::string channel_;
    int torus_;
    bool requireTopology_ = false;

    FiducialCuts* FC_ = nullptr;
    TTree* tree_ = nullptr;

    // UNDER CONSTRUCTION: //

    EventVars ev_;
    RecVars part_;
    RecVars e_;
    RecVars p_;
    RecVars g_;
    GenVars gen_;
    DISVars dis_;
    EPPI0Vars eppi0_;

    struct BranchInfo {
    std::string name;          // e.g. "ev"
    void* address;             // pointer to data
    std::string classType;     // e.g. "EventVars" for TObject-derived structs; empty for primitives
    std::string leafType;      // e.g. "I" or "D" for primitives; empty for TObject structs
};

    std::vector<BranchInfo> masterBranchList_;
    std::unordered_set<std::string> enabledBranches_;

    template <typename T>
    void registerBranchInfo(std::vector<BranchInfo>& list, const std::string& name, T* address, const std::string& classType, const std::string& leafType) {
        list.push_back({name, static_cast<void*>(address), classType, leafType});
    };
    void registerKnownBranchInfo() {

        // TObject-derived structs: provide className, leave leafList empty
        // masterBranchList_.push_back({"ev", ev_, "EventVars", ""});
        // masterBranchList_.push_back({"part", part_, "RecVars", ""});
        // masterBranchList_.push_back({"eppi0", eppi0_, "EPPI0Vars", ""});
        // masterBranchList_.push_back({"dis", dis_, "DISVars", ""});

        // Event info vars
        // registerBranchInfo(masterBranchList_, "runNum",        &ev_.runNum,           "I");
        // registerBranchInfo(masterBranchList_, "eventNum",      &ev_.eventNum,         "I");
        // registerBranchInfo(masterBranchList_, "helicity",      &ev_.helicity,         "I");

        // Reconstructed event vars (inclusive):
        //tree_->Branch("part", &part_);
        // registerBranchInfo(masterBranchList_, "pid",           &part_.pid,             "I");
        // registerBranchInfo(masterBranchList_, "charge",        &part_.charge,          "I");
        // registerBranchInfo(masterBranchList_, "status",        &part_.status,          "I");
        // registerBranchInfo(masterBranchList_, "det",           &part_.det,             "I");
        // registerBranchInfo(masterBranchList_, "sector",        &part_.sector,          "I");
        // registerBranchInfo(masterBranchList_, "p",             &part_.p,               "D");
        // registerBranchInfo(masterBranchList_, "beta",          &part_.beta,            "D");
        // registerBranchInfo(masterBranchList_, "theta",         &part_.theta,           "D");
        // registerBranchInfo(masterBranchList_, "phi",           &part_.phi,             "D");
        // registerBranchInfo(masterBranchList_, "px",            &part_.px,              "D");
        // registerBranchInfo(masterBranchList_, "py",            &part_.py,              "D");
        // registerBranchInfo(masterBranchList_, "pz",            &part_.pz,              "D");
        // registerBranchInfo(masterBranchList_, "vx",            &part_.vx,              "D");
        // registerBranchInfo(masterBranchList_, "vy",            &part_.vy,              "D");
        // registerBranchInfo(masterBranchList_, "vz",            &part_.vz,              "D");
        // registerBranchInfo(masterBranchList_, "chi2pid",       &part_.chi2pid,         "D");
        // registerBranchInfo(masterBranchList_, "time",          &part_.time,            "D");
        // registerBranchInfo(masterBranchList_, "xFT",           &part_.xFT,             "D");
        // registerBranchInfo(masterBranchList_, "yFT",           &part_.yFT,             "D");
        // registerBranchInfo(masterBranchList_, "xDC1",          &part_.xDC1,            "D");
        // registerBranchInfo(masterBranchList_, "yDC1",          &part_.yDC1,            "D");
        // registerBranchInfo(masterBranchList_, "xDC2",          &part_.xDC2,            "D");
        // registerBranchInfo(masterBranchList_, "yDC2",          &part_.yDC2,            "D");
        // registerBranchInfo(masterBranchList_, "xDC3",          &part_.xDC3,            "D");
        // registerBranchInfo(masterBranchList_, "yDC3",          &part_.yDC3,            "D");
        // registerBranchInfo(masterBranchList_, "xPCAL",         &part_.xPCAL,           "D");
        // registerBranchInfo(masterBranchList_, "yPCAL",         &part_.yPCAL,           "D");
        // registerBranchInfo(masterBranchList_, "uPCAL",         &part_.uPCAL,           "D");
        // registerBranchInfo(masterBranchList_, "vPCAL",         &part_.vPCAL,           "D");
        // registerBranchInfo(masterBranchList_, "wPCAL",         &part_.wPCAL,           "D");
        // registerBranchInfo(masterBranchList_, "uECIN",         &part_.uECIN,           "D");
        // registerBranchInfo(masterBranchList_, "vECIN",         &part_.vECIN,           "D");
        // registerBranchInfo(masterBranchList_, "wECIN",         &part_.wECIN,           "D");
        // registerBranchInfo(masterBranchList_, "uECOUT",        &part_.uECOUT,          "D");
        // registerBranchInfo(masterBranchList_, "vECOUT",        &part_.vECOUT,          "D");
        // registerBranchInfo(masterBranchList_, "wECOUT",        &part_.wECOUT,          "D");
        // registerBranchInfo(masterBranchList_, "E_PCAL",        &part_.E_PCAL,          "D");
        // registerBranchInfo(masterBranchList_, "E_ECIN",        &part_.E_ECIN,          "D");
        // registerBranchInfo(masterBranchList_, "E_ECOUT",       &part_.E_ECOUT,         "D");
        // registerBranchInfo(masterBranchList_, "edge_cvt1",     &part_.edge_cvt1,       "D");
        // registerBranchInfo(masterBranchList_, "edge_cvt3",     &part_.edge_cvt3,       "D");
        // registerBranchInfo(masterBranchList_, "edge_cvt5",     &part_.edge_cvt5,       "D");
        // registerBranchInfo(masterBranchList_, "edge_cvt7",     &part_.edge_cvt7,       "D");
        // registerBranchInfo(masterBranchList_, "edge_cvt12",    &part_.edge_cvt12,      "D");
        // registerBranchInfo(masterBranchList_, "theta_cvt",     &part_.theta_cvt,       "D");
        // registerBranchInfo(masterBranchList_, "phi_cvt",       &part_.phi_cvt,         "D");


        // Generated event vars:

        // registerBranchInfo(masterBranchList_, "gen_pid",       &gen_.pid,              "I");
        // registerBranchInfo(masterBranchList_, "gen_p",         &gen_.p,                "D");
        // registerBranchInfo(masterBranchList_, "gen_theta",     &gen_.theta,            "D");
        // registerBranchInfo(masterBranchList_, "gen_phi",       &gen_.phi,              "D");
        // registerBranchInfo(masterBranchList_, "gen_Q2",        &gen_.Q2,               "D");
        // registerBranchInfo(masterBranchList_, "gen_nu",        &gen_.nu,               "D");
        // registerBranchInfo(masterBranchList_, "gen_Xb",        &gen_.Xb,               "D");
        // registerBranchInfo(masterBranchList_, "gen_y",         &gen_.y,                "D");
        // registerBranchInfo(masterBranchList_, "gen_W",         &gen_.W,                "D");
        // registerBranchInfo(masterBranchList_, "gen_t",         &gen_.t,                "D");

        // Electron particle vars:
        // registerBranchInfo(masterBranchList_, "e_p",           &e_.p,                 "D");
        // registerBranchInfo(masterBranchList_, "e_beta",        &e_.beta,              "D");
        // registerBranchInfo(masterBranchList_, "e_theta",       &e_.theta,             "D");
        // registerBranchInfo(masterBranchList_, "e_phi",         &e_.phi,               "D");
        // registerBranchInfo(masterBranchList_, "e_vz",          &e_.vz,                "D");
        // registerBranchInfo(masterBranchList_, "e_xFT",         &e_.xFT,               "D");
        // registerBranchInfo(masterBranchList_, "e_yFT",         &e_.yFT,               "D");
        // registerBranchInfo(masterBranchList_, "e_xDC1",        &e_.xDC1,              "D");
        // registerBranchInfo(masterBranchList_, "e_yDC1",        &e_.yDC1,              "D");
        // registerBranchInfo(masterBranchList_, "e_xDC2",        &e_.xDC2,              "D");
        // registerBranchInfo(masterBranchList_, "e_yDC2",        &e_.yDC2,              "D");
        // registerBranchInfo(masterBranchList_, "e_xDC3",        &e_.xDC3,              "D");
        // registerBranchInfo(masterBranchList_, "e_yDC3",        &e_.yDC3,              "D");
        // registerBranchInfo(masterBranchList_, "e_xPCAL",       &e_.xPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "e_yPCAL",       &e_.yPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "e_uPCAL",       &e_.uPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "e_vPCAL",       &e_.vPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "e_wPCAL",       &e_.wPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "e_uECIN",       &e_.uECIN,             "D");
        // registerBranchInfo(masterBranchList_, "e_vECIN",       &e_.vECIN,             "D");
        // registerBranchInfo(masterBranchList_, "e_wECIN",       &e_.wECIN,             "D");
        // registerBranchInfo(masterBranchList_, "e_uECOUT",      &e_.uECOUT,            "D");
        // registerBranchInfo(masterBranchList_, "e_vECOUT",      &e_.vECOUT,            "D");
        // registerBranchInfo(masterBranchList_, "e_wECOUT",      &e_.wECOUT,            "D");
        // registerBranchInfo(masterBranchList_, "e_E_PCAL",      &e_.E_PCAL,            "D");
        // registerBranchInfo(masterBranchList_, "e_E_ECIN",      &e_.E_ECIN,            "D");
        // registerBranchInfo(masterBranchList_, "e_E_ECOUT",     &e_.E_ECOUT,           "D");
        // registerBranchInfo(masterBranchList_, "e_chi2pid",     &e_.chi2pid,           "D");
        // registerBranchInfo(masterBranchList_, "e_status",      &e_.status,            "I");
        // registerBranchInfo(masterBranchList_, "detEle",        &e_.det,               "I");
        // registerBranchInfo(masterBranchList_, "e_sector",      &e_.sector,            "I");

        // // Proton particle vars:
        // registerBranchInfo(masterBranchList_, "p_p",           &p_.p,                 "D");
        // registerBranchInfo(masterBranchList_, "p_beta",        &p_.beta,              "D");
        // registerBranchInfo(masterBranchList_, "p_theta",       &p_.theta,             "D");
        // registerBranchInfo(masterBranchList_, "p_phi",         &p_.phi,               "D");
        // registerBranchInfo(masterBranchList_, "p_vz",          &p_.vz,                "D");
        // registerBranchInfo(masterBranchList_, "p_xDC1",        &p_.xDC1,              "D");
        // registerBranchInfo(masterBranchList_, "p_yDC1",        &p_.yDC1,              "D");
        // registerBranchInfo(masterBranchList_, "p_xDC2",        &p_.xDC2,              "D");
        // registerBranchInfo(masterBranchList_, "p_yDC2",        &p_.yDC2,              "D");
        // registerBranchInfo(masterBranchList_, "p_xDC3",        &p_.xDC3,              "D");
        // registerBranchInfo(masterBranchList_, "p_yDC3",        &p_.yDC3,              "D");
        // registerBranchInfo(masterBranchList_, "p_chi2pid",     &p_.chi2pid,           "D");
        // registerBranchInfo(masterBranchList_, "p_edge_cvt1",   &p_.edge_cvt1,         "D");
        // registerBranchInfo(masterBranchList_, "p_edge_cvt3",   &p_.edge_cvt3,         "D");
        // registerBranchInfo(masterBranchList_, "p_edge_cvt5",   &p_.edge_cvt5,         "D");
        // registerBranchInfo(masterBranchList_, "p_edge_cvt7",   &p_.edge_cvt7,         "D");
        // registerBranchInfo(masterBranchList_, "p_edge_cvt12",  &p_.edge_cvt12,        "D");
        // registerBranchInfo(masterBranchList_, "p_theta_cvt",   &p_.theta_cvt,         "D");
        // registerBranchInfo(masterBranchList_, "p_phi_cvt",     &p_.phi_cvt,           "D");
        // registerBranchInfo(masterBranchList_, "p_status",      &p_.status,            "I");
        // registerBranchInfo(masterBranchList_, "detPro",        &p_.det,               "I");
        // registerBranchInfo(masterBranchList_, "p_sector",      &p_.sector,            "I");

        // // Photon particle vars:
        // registerBranchInfo(masterBranchList_, "g_p",           &g_.p,                 "D");
        // registerBranchInfo(masterBranchList_, "g_theta",       &g_.theta,             "D");
        // registerBranchInfo(masterBranchList_, "g_phi",         &g_.phi,               "D");
        // registerBranchInfo(masterBranchList_, "g_vz",          &g_.vz,                "D");
        // registerBranchInfo(masterBranchList_, "g_xFT",         &g_.xFT,               "D");
        // registerBranchInfo(masterBranchList_, "g_yFT",         &g_.yFT,               "D");
        // registerBranchInfo(masterBranchList_, "g_xPCAL",       &g_.xPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "g_yPCAL",       &g_.yPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "g_uPCAL",       &g_.uPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "g_vPCAL",       &g_.vPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "g_wPCAL",       &g_.wPCAL,             "D");
        // registerBranchInfo(masterBranchList_, "g_uECIN",       &g_.uECIN,             "D");
        // registerBranchInfo(masterBranchList_, "g_vECIN",       &g_.vECIN,             "D");
        // registerBranchInfo(masterBranchList_, "g_wECIN",       &g_.wECIN,             "D");
        // registerBranchInfo(masterBranchList_, "g_uECOUT",      &g_.uECOUT,            "D");
        // registerBranchInfo(masterBranchList_, "g_vECOUT",      &g_.vECOUT,            "D");
        // registerBranchInfo(masterBranchList_, "g_wECOUT",      &g_.wECOUT,            "D");
        // registerBranchInfo(masterBranchList_, "g_E_PCAL",      &g_.E_PCAL,            "D");
        // registerBranchInfo(masterBranchList_, "g_E_ECIN",      &g_.E_ECIN,            "D");
        // registerBranchInfo(masterBranchList_, "g_E_ECOUT",     &g_.E_ECOUT,           "D");
        // registerBranchInfo(masterBranchList_, "g_chi2pid",     &g_.chi2pid,           "D");
        // registerBranchInfo(masterBranchList_, "g_status",      &g_.status,            "I");
        // registerBranchInfo(masterBranchList_, "detPho",        &g_.det,               "I");
        // registerBranchInfo(masterBranchList_, "g_sector",      &g_.sector,            "I");

        // Pi0 & exclusivity vars:
        // registerBranchInfo(masterBranchList_, "m_gg",          &eppi0_.m_gg,          "D");
        // registerBranchInfo(masterBranchList_, "pi0_p",         &eppi0_.pi0_p,         "D");
        // registerBranchInfo(masterBranchList_, "pi0_theta",     &eppi0_.pi0_theta,     "D");
        // registerBranchInfo(masterBranchList_, "pi0_phi",       &eppi0_.pi0_phi,       "D");
        // registerBranchInfo(masterBranchList_, "t",             &eppi0_.t,             "D");
        // registerBranchInfo(masterBranchList_, "m2_miss",       &eppi0_.m2_miss,       "D");
        // registerBranchInfo(masterBranchList_, "m2_epX",        &eppi0_.m2_epX,        "D");
        // registerBranchInfo(masterBranchList_, "m2_epi0X",      &eppi0_.m2_epi0X,      "D");
        // registerBranchInfo(masterBranchList_, "E_miss",        &eppi0_.E_miss,        "D");
        // registerBranchInfo(masterBranchList_, "px_miss",       &eppi0_.px_miss,       "D");
        // registerBranchInfo(masterBranchList_, "py_miss",       &eppi0_.py_miss,       "D");
        // registerBranchInfo(masterBranchList_, "pz_miss",       &eppi0_.pz_miss,       "D");
        // registerBranchInfo(masterBranchList_, "pT_miss",       &eppi0_.pT_miss,       "D");
        // registerBranchInfo(masterBranchList_, "trentoPhi",     &eppi0_.trentoPhi,     "D");
        // registerBranchInfo(masterBranchList_, "detPi0",        &eppi0_.detPi0,        "I");
        // registerBranchInfo(masterBranchList_, "numElectrons",  &eppi0_.numElectrons,  "I");
        // registerBranchInfo(masterBranchList_, "numProtons",    &eppi0_.numProtons,    "I");
        // registerBranchInfo(masterBranchList_, "numPhotons",    &eppi0_.numPhotons,    "I");

        // DIS vars
        // registerBranchInfo(masterBranchList_, "Q2",            &dis_.Q2,               "D");
        // registerBranchInfo(masterBranchList_, "nu",            &dis_.nu,               "D");
        // registerBranchInfo(masterBranchList_, "Xb",            &dis_.Xb,               "D");
        // registerBranchInfo(masterBranchList_, "y",             &dis_.y,                "D");
        // registerBranchInfo(masterBranchList_, "W",             &dis_.W,                "D");

    };
    void attachEnabledBranches(TTree* tree) {
        for (const auto& brInfo : masterBranchList_) {
            if (enabledBranches_.count(brInfo.name)) {
                if (!brInfo.classType.empty()) {
                    // TObject-derived struct
                    tree->Branch(brInfo.name.c_str(), brInfo.classType.c_str(), brInfo.address);
                } else if (!brInfo.leafType.empty()) {
                    // Primitive data type
                    std::string branchDef = brInfo.name + "/" + brInfo.leafType;
                    tree->Branch(brInfo.name.c_str(), brInfo.address, branchDef.c_str());
                } else {
                    // Unexpected case: no classType and no leafType
                    // Handle error or warning as needed
                    std::cerr << "Warning: Branch Info " << brInfo.name << " has no className or leafList." << std::endl;
                }
            }
        }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Variables to be stored in tree are declared here, linked to branches in rootTree(), and filled/updated in processEvent(): ///
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Inclusive:
    double   chi2pid_;
    int      pid_, charge_;

    // MC vars:
    double   e_pgen_, p_pgen_, g_pgen_;
    double   Q2_gen_, nu_gen_, Xb_gen_, y_gen_, W_gen_, t_gen_;

    // Electron:
    double   e_p_, e_beta_, e_theta_, e_phi_, e_vz_;
    double   e_xFT_, e_yFT_;
    double   e_xDC1_, e_yDC1_, e_xDC2_, e_yDC2_, e_xDC3_, e_yDC3_;
    double   e_xPCAL_, e_yPCAL_, e_uPCAL_, e_vPCAL_, e_wPCAL_, e_uECIN_, e_vECIN_, e_wECIN_, e_uECOUT_, e_vECOUT_, e_wECOUT_;
    double   e_E_PCAL_, e_E_ECIN_, e_E_ECOUT_;
    double   e_chi2pid_;
    int      e_status_, detEle_, e_region_, e_sector_, nphe_;
    
    // Proton:
    double   p_p_, p_beta_, p_theta_, p_phi_, p_vz_;
    double   p_xDC1_, p_yDC1_, p_xDC2_, p_yDC2_, p_xDC3_, p_yDC3_;
    double   p_chi2pid_;
    int      p_status_, detPro_, p_sector_;
    double   p_edge_cvt1_, p_edge_cvt3_, p_edge_cvt5_, p_edge_cvt7_, p_edge_cvt12_;
    double   p_theta_cvt_, p_phi_cvt_;

    // Photon(s):
    double   g_p_, g_theta_, g_phi_, g_vz_;
    double   g_xPCAL_, g_yPCAL_, g_uPCAL_, g_vPCAL_, g_wPCAL_, g_uECIN_, g_vECIN_, g_wECIN_, g_uECOUT_, g_vECOUT_, g_wECOUT_;
    double   g_E_PCAL_, g_E_ECIN_, g_E_ECOUT_;
    double   g_chi2pid_;
    int      g_status_, detPho_, g_sector_;

    // Pi0:
    double   m_gg_, pi0_p_, pi0_theta_, pi0_phi_; // pi0 quantities gotten from photon pair
    int      detPi0_;

    // Event information to store:
    int      runNum_, eventNum_, helicity_;
    int      numElectrons_, numProtons_, numPhotons_;
    
    // DIS quantities:
    double   Q2_, nu_, Xb_, y_, W_;

    // Exclusivity & derived quantities:
    double   t_;
    double   m2_miss_, m2_epX_, m2_epi0X_,E_miss_, E2_miss_;
    double   px_miss_, py_miss_, pz_miss_, pT_miss_;
    double   trentoPhi_, deltaPhi_;

};

