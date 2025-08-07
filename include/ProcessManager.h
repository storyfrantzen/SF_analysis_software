#pragma once

#include "nlohmann/json.hpp"
#include "FiducialCuts.h"  // Note: ProcessManager has its own FC object during filtering
#include "Kinematics.h" // Note: ProcessManager will use kinematics class for computing derived quantities
#include "PhysicalConstants.h" // contains useful constants
#include "clas12reader.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <string>
#include <cstdint>  // for int8_t
#include <iomanip>  // for std::put_time
#include <sstream>  // for std::ostringstream
#include <ctime>    // for std::time_t, std::localtime

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
    struct GenVars {;};
    struct RecVars {
        double p, beta, theta, phi, vz, chi2pid; 
        int pid, charge, status, det, sector;

        void fill(clas12::region_particle* rec) {
            p = rec->getP();
            beta = rec->par()->getBeta();
            theta = rec->getTheta();
            phi = rec->getPhi();
            vz = rec->par()->getVz();
            chi2pid = rec->par()->getChi2Pid();

            pid = rec->getPid();
            charge = rec->par()->getCharge();
            status = rec->par()->getStatus();
            det = ProcessManager::getDetector(status);
            sector = rec->getSector();
        }
    };
    struct ElectronVars { double p = 0, theta = 0, phi = 0; };
    struct ProtonVars   { double p = 0; };
    struct PhotonVars   { double px = 0, py = 0, pz = 0; };

    struct BranchInfo {
        std::string name;
        void* address;
        std::string type; // e.g., "D" for double
    };

    template <typename T>
    void registerBranchInfo(std::vector<BranchInfo>& list, const std::string& name, T* address, const std::string& type) {
        list.push_back({name, static_cast<void*>(address), type});
    };
    void registerAllBranchInfo() {
        // Inclusive
        registerBranchInfo(masterBranchList_, "p", &rec_.p, "D");
        registerBranchInfo(masterBranchList_, "beta", &rec_.beta, "D");
        registerBranchInfo(masterBranchList_, "theta", &rec_.theta, "D");
        registerBranchInfo(masterBranchList_, "phi", &rec_.phi, "D");
        registerBranchInfo(masterBranchList_, "vz", &rec_.vz, "D");
        registerBranchInfo(masterBranchList_, "chi2pid", &rec_.chi2pid, "D");
        registerBranchInfo(masterBranchList_, "pid", &rec_.pid, "I");
        registerBranchInfo(masterBranchList_, "charge", &rec_.charge, "I");
        registerBranchInfo(masterBranchList_, "status", &rec_.status, "I");
        registerBranchInfo(masterBranchList_, "sector", &rec_.sector, "I");

        // MC vars
        registerBranchInfo(masterBranchList_, "e_pgen", &e_pgen_, "D");
        registerBranchInfo(masterBranchList_, "p_pgen", &p_pgen_, "D");
        registerBranchInfo(masterBranchList_, "g_pgen", &g_pgen_, "D");
        registerBranchInfo(masterBranchList_, "Q2_gen", &Q2_gen_, "D");
        registerBranchInfo(masterBranchList_, "nu_gen", &nu_gen_, "D");
        registerBranchInfo(masterBranchList_, "Xb_gen", &Xb_gen_, "D");
        registerBranchInfo(masterBranchList_, "y_gen", &y_gen_, "D");
        registerBranchInfo(masterBranchList_, "W_gen", &W_gen_, "D");
        registerBranchInfo(masterBranchList_, "t_gen", &t_gen_, "D");

        // Electron
        registerBranchInfo(masterBranchList_, "e_p", &e_p_, "D");
        registerBranchInfo(masterBranchList_, "e_beta", &e_beta_, "D");
        registerBranchInfo(masterBranchList_, "e_theta", &e_theta_, "D");
        registerBranchInfo(masterBranchList_, "e_phi", &e_phi_, "D");
        registerBranchInfo(masterBranchList_, "e_vz", &e_vz_, "D");
        registerBranchInfo(masterBranchList_, "e_xFT", &e_xFT_, "D");
        registerBranchInfo(masterBranchList_, "e_yFT", &e_yFT_, "D");
        registerBranchInfo(masterBranchList_, "e_xDC1", &e_xDC1_, "D");
        registerBranchInfo(masterBranchList_, "e_yDC1", &e_yDC1_, "D");
        registerBranchInfo(masterBranchList_, "e_xDC2", &e_xDC2_, "D");
        registerBranchInfo(masterBranchList_, "e_yDC2", &e_yDC2_, "D");
        registerBranchInfo(masterBranchList_, "e_xDC3", &e_xDC3_, "D");
        registerBranchInfo(masterBranchList_, "e_yDC3", &e_yDC3_, "D");
        registerBranchInfo(masterBranchList_, "e_xPCAL", &e_xPCAL_, "D");
        registerBranchInfo(masterBranchList_, "e_yPCAL", &e_yPCAL_, "D");
        registerBranchInfo(masterBranchList_, "e_uPCAL", &e_uPCAL_, "D");
        registerBranchInfo(masterBranchList_, "e_vPCAL", &e_vPCAL_, "D");
        registerBranchInfo(masterBranchList_, "e_wPCAL", &e_wPCAL_, "D");
        registerBranchInfo(masterBranchList_, "e_uECIN", &e_uECIN_, "D");
        registerBranchInfo(masterBranchList_, "e_vECIN", &e_vECIN_, "D");
        registerBranchInfo(masterBranchList_, "e_wECIN", &e_wECIN_, "D");
        registerBranchInfo(masterBranchList_, "e_uECOUT", &e_uECOUT_, "D");
        registerBranchInfo(masterBranchList_, "e_vECOUT", &e_vECOUT_, "D");
        registerBranchInfo(masterBranchList_, "e_wECOUT", &e_wECOUT_, "D");
        registerBranchInfo(masterBranchList_, "e_E_PCAL", &e_E_PCAL_, "D");
        registerBranchInfo(masterBranchList_, "e_E_ECIN", &e_E_ECIN_, "D");
        registerBranchInfo(masterBranchList_, "e_E_ECOUT", &e_E_ECOUT_, "D");
        registerBranchInfo(masterBranchList_, "e_chi2pid", &e_chi2pid_, "D");
        registerBranchInfo(masterBranchList_, "e_status", &e_status_, "I");
        registerBranchInfo(masterBranchList_, "detEle", &detEle_, "I");
        registerBranchInfo(masterBranchList_, "e_region", &e_region_, "I");
        registerBranchInfo(masterBranchList_, "e_sector", &e_sector_, "I");
        registerBranchInfo(masterBranchList_, "nphe", &nphe_, "I");

        // Proton
        registerBranchInfo(masterBranchList_, "p_p", &p_p_, "D");
        registerBranchInfo(masterBranchList_, "p_beta", &p_beta_, "D");
        registerBranchInfo(masterBranchList_, "p_theta", &p_theta_, "D");
        registerBranchInfo(masterBranchList_, "p_phi", &p_phi_, "D");
        registerBranchInfo(masterBranchList_, "p_vz", &p_vz_, "D");
        registerBranchInfo(masterBranchList_, "p_xDC1", &p_xDC1_, "D");
        registerBranchInfo(masterBranchList_, "p_yDC1", &p_yDC1_, "D");
        registerBranchInfo(masterBranchList_, "p_xDC2", &p_xDC2_, "D");
        registerBranchInfo(masterBranchList_, "p_yDC2", &p_yDC2_, "D");
        registerBranchInfo(masterBranchList_, "p_xDC3", &p_xDC3_, "D");
        registerBranchInfo(masterBranchList_, "p_yDC3", &p_yDC3_, "D");
        registerBranchInfo(masterBranchList_, "p_chi2pid", &p_chi2pid_, "D");
        registerBranchInfo(masterBranchList_, "p_status", &p_status_, "I");
        registerBranchInfo(masterBranchList_, "detPro", &detPro_, "I");
        registerBranchInfo(masterBranchList_, "p_sector", &p_sector_, "I");
        registerBranchInfo(masterBranchList_, "p_edge_cvt1", &p_edge_cvt1_, "D");
        registerBranchInfo(masterBranchList_, "p_edge_cvt3", &p_edge_cvt3_, "D");
        registerBranchInfo(masterBranchList_, "p_edge_cvt5", &p_edge_cvt5_, "D");
        registerBranchInfo(masterBranchList_, "p_edge_cvt7", &p_edge_cvt7_, "D");
        registerBranchInfo(masterBranchList_, "p_edge_cvt12", &p_edge_cvt12_, "D");
        registerBranchInfo(masterBranchList_, "p_theta_cvt", &p_theta_cvt_, "D");
        registerBranchInfo(masterBranchList_, "p_phi_cvt", &p_phi_cvt_, "D");

        // Photon
        registerBranchInfo(masterBranchList_, "g_p", &g_p_, "D");
        registerBranchInfo(masterBranchList_, "g_theta", &g_theta_, "D");
        registerBranchInfo(masterBranchList_, "g_phi", &g_phi_, "D");
        registerBranchInfo(masterBranchList_, "g_vz", &g_vz_, "D");
        registerBranchInfo(masterBranchList_, "g_xPCAL", &g_xPCAL_, "D");
        registerBranchInfo(masterBranchList_, "g_yPCAL", &g_yPCAL_, "D");
        registerBranchInfo(masterBranchList_, "g_uPCAL", &g_uPCAL_, "D");
        registerBranchInfo(masterBranchList_, "g_vPCAL", &g_vPCAL_, "D");
        registerBranchInfo(masterBranchList_, "g_wPCAL", &g_wPCAL_, "D");
        registerBranchInfo(masterBranchList_, "g_uECIN", &g_uECIN_, "D");
        registerBranchInfo(masterBranchList_, "g_vECIN", &g_vECIN_, "D");
        registerBranchInfo(masterBranchList_, "g_wECIN", &g_wECIN_, "D");
        registerBranchInfo(masterBranchList_, "g_uECOUT", &g_uECOUT_, "D");
        registerBranchInfo(masterBranchList_, "g_vECOUT", &g_vECOUT_, "D");
        registerBranchInfo(masterBranchList_, "g_wECOUT", &g_wECOUT_, "D");
        registerBranchInfo(masterBranchList_, "g_E_PCAL", &g_E_PCAL_, "D");
        registerBranchInfo(masterBranchList_, "g_E_ECIN", &g_E_ECIN_, "D");
        registerBranchInfo(masterBranchList_, "g_E_ECOUT", &g_E_ECOUT_, "D");
        registerBranchInfo(masterBranchList_, "g_chi2pid", &g_chi2pid_, "D");
        registerBranchInfo(masterBranchList_, "g_status", &g_status_, "I");
        registerBranchInfo(masterBranchList_, "detPho", &detPho_, "I");
        registerBranchInfo(masterBranchList_, "g_sector", &g_sector_, "I");

        // Pi0
        registerBranchInfo(masterBranchList_, "m_gg", &m_gg_, "D");
        registerBranchInfo(masterBranchList_, "pi0_p", &pi0_p_, "D");
        registerBranchInfo(masterBranchList_, "pi0_theta", &pi0_theta_, "D");
        registerBranchInfo(masterBranchList_, "pi0_phi", &pi0_phi_, "D");
        registerBranchInfo(masterBranchList_, "detPi0", &detPi0_, "I");

        // Pi+
        registerBranchInfo(masterBranchList_, "pip_p", &pip_p_, "D");
        registerBranchInfo(masterBranchList_, "pip_theta", &pip_theta_, "D");
        registerBranchInfo(masterBranchList_, "pip_phi", &pip_phi_, "D");
        registerBranchInfo(masterBranchList_, "pip_xDC1", &pip_xDC1_, "D");
        registerBranchInfo(masterBranchList_, "pip_yDC1", &pip_yDC1_, "D");
        registerBranchInfo(masterBranchList_, "pip_xDC2", &pip_xDC2_, "D");
        registerBranchInfo(masterBranchList_, "pip_yDC2", &pip_yDC2_, "D");
        registerBranchInfo(masterBranchList_, "pip_xDC3", &pip_xDC3_, "D");
        registerBranchInfo(masterBranchList_, "pip_yDC3", &pip_yDC3_, "D");
        registerBranchInfo(masterBranchList_, "pip_E_PCAL", &pip_E_PCAL_, "D");
        registerBranchInfo(masterBranchList_, "pip_E_ECIN", &pip_E_ECIN_, "D");
        registerBranchInfo(masterBranchList_, "pip_E_ECOUT", &pip_E_ECOUT_, "D");
        registerBranchInfo(masterBranchList_, "detPip", &detPip_, "I");
        registerBranchInfo(masterBranchList_, "pip_sector", &pip_sector_, "I");

        // Pi-
        registerBranchInfo(masterBranchList_, "pim_p", &pim_p_, "D");
        registerBranchInfo(masterBranchList_, "pim_theta", &pim_theta_, "D");
        registerBranchInfo(masterBranchList_, "pim_phi", &pim_phi_, "D");
        registerBranchInfo(masterBranchList_, "pim_xDC1", &pim_xDC1_, "D");
        registerBranchInfo(masterBranchList_, "pim_yDC1", &pim_yDC1_, "D");
        registerBranchInfo(masterBranchList_, "pim_xDC2", &pim_xDC2_, "D");
        registerBranchInfo(masterBranchList_, "pim_yDC2", &pim_yDC2_, "D");
        registerBranchInfo(masterBranchList_, "pim_xDC3", &pim_xDC3_, "D");
        registerBranchInfo(masterBranchList_, "pim_yDC3", &pim_yDC3_, "D");
        registerBranchInfo(masterBranchList_, "pim_E_PCAL", &pim_E_PCAL_, "D");
        registerBranchInfo(masterBranchList_, "pim_E_ECIN", &pim_E_ECIN_, "D");
        registerBranchInfo(masterBranchList_, "pim_E_ECOUT", &pim_E_ECOUT_, "D");
        registerBranchInfo(masterBranchList_, "detPim", &detPim_, "I");
        registerBranchInfo(masterBranchList_, "pim_sector", &pim_sector_, "I");

        // Event info
        registerBranchInfo(masterBranchList_, "runNum", &runNum_, "I");
        registerBranchInfo(masterBranchList_, "eventNum", &eventNum_, "I");
        registerBranchInfo(masterBranchList_, "helicity", &helicity_, "I");
        registerBranchInfo(masterBranchList_, "numElectrons", &numElectrons_, "I");
        registerBranchInfo(masterBranchList_, "numProtons", &numProtons_, "I");
        registerBranchInfo(masterBranchList_, "numPhotons", &numPhotons_, "I");
        registerBranchInfo(masterBranchList_, "numPip", &numPip_, "I");
        registerBranchInfo(masterBranchList_, "numPim", &numPim_, "I");

        // DIS
        registerBranchInfo(masterBranchList_, "Q2", &Q2_, "D");
        registerBranchInfo(masterBranchList_, "nu", &nu_, "D");
        registerBranchInfo(masterBranchList_, "Xb", &Xb_, "D");
        registerBranchInfo(masterBranchList_, "y", &y_, "D");
        registerBranchInfo(masterBranchList_, "W", &W_, "D");

        // Exclusivity
        registerBranchInfo(masterBranchList_, "t", &t_, "D");
        registerBranchInfo(masterBranchList_, "m2_miss", &m2_miss_, "D");
        registerBranchInfo(masterBranchList_, "m2_epX", &m2_epX_, "D");
        registerBranchInfo(masterBranchList_, "m2_epi0X", &m2_epi0X_, "D");
        registerBranchInfo(masterBranchList_, "E_miss", &E_miss_, "D");
        registerBranchInfo(masterBranchList_, "E2_miss", &E2_miss_, "D");
        registerBranchInfo(masterBranchList_, "px_miss", &px_miss_, "D");
        registerBranchInfo(masterBranchList_, "py_miss", &py_miss_, "D");
        registerBranchInfo(masterBranchList_, "pz_miss", &pz_miss_, "D");
        registerBranchInfo(masterBranchList_, "pT_miss", &pT_miss_, "D");
        registerBranchInfo(masterBranchList_, "trentoPhi", &trentoPhi_, "D");
        registerBranchInfo(masterBranchList_, "deltaPhi", &deltaPhi_, "D");

    };
    void attachEnabledBranches(TTree* tree) {
        for (const auto& b : masterBranchList_) {
            if (enabledBranches_.count(b.name)) {
                tree->Branch(b.name.c_str(), b.address, (b.name + "/" + b.type).c_str());
            }
        }
    };
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Variables to be stored in tree are declared here, linked to branches in rootTree(), and filled/updated in processEvent(): ///
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    RecVars rec_;
    ElectronVars e_;
    ProtonVars p_;
    PhotonVars g_;

    std::vector<BranchInfo> masterBranchList_;
    std::unordered_set<std::string> enabledBranches_;

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

    // Pi +/-:
    double   pip_p_, pip_theta_, pip_phi_;
    double   pip_xDC1_, pip_yDC1_, pip_xDC2_, pip_yDC2_, pip_xDC3_, pip_yDC3_;
    double   pip_E_PCAL_, pip_E_ECIN_, pip_E_ECOUT_;
    int      detPip_, pip_sector_;

    double   pim_p_, pim_theta_, pim_phi_;
    double   pim_xDC1_, pim_yDC1_, pim_xDC2_, pim_yDC2_, pim_xDC3_, pim_yDC3_;
    double   pim_E_PCAL_, pim_E_ECIN_, pim_E_ECOUT_;
    int      detPim_, pim_sector_;

    // Event information to store:
    int      runNum_, eventNum_, helicity_;
    int      numElectrons_, numProtons_, numPhotons_, numPip_, numPim_;
    
    // DIS quantities:
    double   Q2_, nu_, Xb_, y_, W_;

    // Exclusivity & derived quantities:
    double   t_;
    double   m2_miss_, m2_epX_, m2_epi0X_,E_miss_, E2_miss_;
    double   px_miss_, py_miss_, pz_miss_, pT_miss_;
    double   trentoPhi_, deltaPhi_;

};

