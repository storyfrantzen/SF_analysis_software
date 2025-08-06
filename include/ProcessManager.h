#pragma once

#include "FiducialCuts.h"  // Note: ProcessManager will have a member copy of the FC during filtering
#include "Kinematics.h" // Note: ProcessManager will use kinematics class for computing derived quantities
#include "clas12reader.h"
#include "TFile.h"
#include "TTree.h"
#include <memory>
#include <string>

class ProcessManager {
public:

    // Setters defined here: //
    void setChannel(const std::string channel) {channel_ = channel;}
    void setTorus(const int torus) {torus_ = torus;}
    void setTopology(const std::vector<std::string>& topology) {
    if (topology.size() != 2) {
        std::cerr << "[Warning] setTopology: Search for 2 elements was unsuccessful." << std::endl;
        requireTopology_ = false;
        return;
    }

    const std::string& pho = topology[0];
    const std::string& pro = topology[1];

    if (pho != "FT" && pho != "FD" && pho != "CD") {
        std::cerr << "[ERROR] setTopology: void topology. Invalid photon detector: " << pho 
                  << " (expected FT, FD, or CD)" << std::endl;
        requireTopology_ = false;
        return;
    }

    if (pro != "FD" && pro != "CD") {
        std::cerr << "[ERROR] setTopology: void topology. Invalid proton detector: " << pro 
                  << " (expected FD or CD)" << std::endl;
        requireTopology_ = false;
        return;
    }

    requireTopology_ = true;

    detPho_ = (pho == "FT") ? 0 : (pho == "FD") ? 1 : 2;
    detPro_ = (pro == "FD") ? 1 : 2;
}
    void setEbeam(const float Ebeam) {Ebeam_ = Ebeam;}
    void setFiducialCuts(FiducialCuts& fiducialCuts) {fiducialCuts_ = &fiducialCuts;}

    // Getters defined here: //
    int eventsProcessed() {return eventsProcessed_;}
    int numFills() {return numFills_;}

    // LITTLE FUNCTIONS: //
    std::string currentTimestamp() const; 
    std::string makeFilename() const;
    bool passesVertexCut(clas12::region_particle* p, const int zmin=-8, const int zmax=2);
    bool passesDiagECALCut(clas12::region_particle* ele);
    int  getDetector(int status);

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
    void processEvent(clas12::clas12reader& c12);
    void finalize(const std::string& output_file);

private:

    int eventsProcessed_ = 0;
    int numFills_ = 0;

    double Ebeam_;
    std::string channel_;
    int torus_;
    bool requireTopology_ = false;

    FiducialCuts* fiducialCuts_ = nullptr;

    TFile* outFile_ = nullptr;
    TTree* tree_ = nullptr;
    
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

