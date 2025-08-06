#include "ProcessManager.h"
#include "Kinematics.h" // Note: ProcessManager will use kinematics class for computing derived quantities
#include "TLorentzVector.h"
#include "PhysicalConstants.h" // contains useful constants
#include <cmath>
#include <cstdint>  // for int8_t
#include <iomanip>  // for std::put_time
#include <sstream>  // for std::ostringstream
#include <ctime>    // for std::time_t, std::localtime

using namespace clas12;

//////////////////////////////////////////////////////////
// Little helper functions called by bigger functions: ///
//////////////////////////////////////////////////////////

std::string ProcessManager::currentTimestamp() const {
    auto now = std::chrono::system_clock::now();
    std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&tt), "%Y%m%d_%H%M%S");
    return ss.str();
}

std::string ProcessManager::makeFilename() const {
    std::stringstream ss;
    ss << Ebeam_ << "_" << torus_ << "_" << channel_ << "_" << eventsProcessed_ << "_" << currentTimestamp() << ".root";
    return ss.str();
}

// Currently in use. Default is -8 <= z <= 2.
bool ProcessManager::passesVertexCut(clas12::region_particle* p, const int zmin, const int zmax) { 
    return p->par()->getVz() >= zmin && p->par()->getVz() <= zmax; 
}

// Currently not in use. RGK analysis note says 
bool ProcessManager::passesDiagECALCut(clas12::region_particle* ele) {
    return true;
}

// Based on particle status, returns int representation of FT (0), FD (1), or CD (2)
int ProcessManager::getDetector(int status) {
    const int absStatus = std::abs(status);

    if (absStatus >= 1000 && absStatus < 2000) return 0; // FT
    if (absStatus >= 2000 && absStatus < 4000) return 1; // FD
    if (absStatus >= 4000 && absStatus < 5000) return 2; // CD

    return -999; // Unknown detector
}

//////////////////////////////////////////////////////////
// Bigger functions: //
//////////////////////////////////////////////////////////

// creates & links electron branches, i.e., with corresponding member variables (e.g. `var_`)
void ProcessManager::writeEleBranches(int info) {
    if (info == 1) {
        tree_->Branch("e_p",          &e_p_,           "e_p/D");
        tree_->Branch("e_beta",       &e_beta_,        "e_beta/D");
        tree_->Branch("e_theta",      &e_theta_,       "e_theta/D");
        tree_->Branch("e_phi",        &e_phi_,         "e_phi/D");
        tree_->Branch("e_vz",         &e_vz_,          "e_vz/D");
        tree_->Branch("e_sector",     &e_sector_,      "e_sector/I");
        tree_->Branch("e_status",     &e_status_,      "e_status/I");
        tree_->Branch("e_chi2pid",    &e_chi2pid_,     "e_chi2pid/D");
        tree_->Branch("detEle",       &detEle_,        "detEle/I");
        // DIS branches:
        tree_->Branch("Q2",           &Q2_,            "Q2/D");
        tree_->Branch("nu",           &nu_,            "nu/D");
        tree_->Branch("Xb",           &Xb_,            "Xb/D");
        tree_->Branch("y",            &y_,             "y/D");
        tree_->Branch("W",            &W_,             "W/D");
        // Misc. branches:
        tree_->Branch("e_xDC1",       &e_xDC1_,        "e_xDC1/D");
        tree_->Branch("e_yDC1",       &e_yDC1_,        "e_yDC1/D");
        tree_->Branch("e_xDC2",       &e_xDC2_,        "e_xDC2/D");
        tree_->Branch("e_yDC2",       &e_yDC2_,        "e_yDC2/D");
        tree_->Branch("e_xDC3",       &e_xDC3_,        "e_xDC3/D");
        tree_->Branch("e_yDC3",       &e_yDC3_,        "e_yDC3/D");
        tree_->Branch("e_uPCAL",      &e_uPCAL_,       "e_uPCAL/D");
        tree_->Branch("e_vPCAL",      &e_vPCAL_,       "e_vPCAL/D");
        tree_->Branch("e_wPCAL",      &e_wPCAL_,       "e_wPCAL/D");
        tree_->Branch("e_uECIN",      &e_uECIN_,       "e_uECIN/D");
        tree_->Branch("e_vECIN",      &e_vECIN_,       "e_vECIN/D");
        tree_->Branch("e_wECIN",      &e_wECIN_,       "e_wECIN/D");
        tree_->Branch("e_uECOUT",     &e_uECOUT_,      "e_uECOUT/D");
        tree_->Branch("e_vECOUT",     &e_vECOUT_,      "e_vECOUT/D");
        tree_->Branch("e_wECOUT",     &e_wECOUT_,      "e_wECOUT/D");
        tree_->Branch("e_E_PCAL",     &e_E_PCAL_,      "e_E_PCAL/D");
        tree_->Branch("e_E_ECIN",     &e_E_ECIN_,      "e_E_ECIN/D");
        tree_->Branch("e_E_ECOUT",    &e_E_ECOUT_,     "e_E_ECOUT/D");
    }
}

// creates & links proton branches
void ProcessManager::writeProBranches(int info) {
    if (info > 0) {
        tree_->Branch("p_p",          &p_p_,           "p_p/D");
        tree_->Branch("p_theta",      &p_theta_,       "p_theta/D");
        tree_->Branch("p_beta",       &p_beta_,        "p_beta/D");
        tree_->Branch("p_phi",        &p_phi_,         "p_phi/D");
        tree_->Branch("p_vz",         &p_vz_,          "p_vz/D");
        tree_->Branch("p_sector",     &p_sector_,      "p_sector/I");
        tree_->Branch("p_status",     &p_status_,      "p_status/I");
        tree_->Branch("p_chi2pid",    &p_chi2pid_,     "p_chi2pid/D");
        tree_->Branch("detPro",       &detPro_,        "detPro/I");
        tree_->Branch("t",            &t_,             "t/D");
        tree_->Branch("p_xDC1",       &p_xDC1_,        "p_xDC1/D");
        tree_->Branch("p_yDC1",       &p_yDC1_,        "p_yDC1/D");
        tree_->Branch("p_xDC2",       &p_xDC2_,        "p_xDC2/D");
        tree_->Branch("p_yDC2",       &p_yDC2_,        "p_yDC2/D");
        tree_->Branch("p_xDC3",       &p_xDC3_,        "p_xDC3/D");
        tree_->Branch("p_yDC3",       &p_yDC3_,        "p_yDC3/D");
        tree_->Branch("p_edge_cvt1",  &p_edge_cvt1_,   "p_edge_cvt1/D");
        tree_->Branch("p_edge_cvt3",  &p_edge_cvt3_,   "p_edge_cvt3/D");
        tree_->Branch("p_edge_cvt5",  &p_edge_cvt5_,   "p_edge_cvt5/D");
        tree_->Branch("p_edge_cvt7",  &p_edge_cvt7_,   "p_edge_cvt7/D");
        tree_->Branch("p_edge_cvt12", &p_edge_cvt12_,  "p_edge_cvt12/D");
        tree_->Branch("p_theta_cvt",  &p_theta_cvt_,   "p_theta_cvt/D");
        tree_->Branch("p_phi_cvt",    &p_phi_cvt_,     "p_phi_cvt/D");
    }
}

// create & links photon branches
void ProcessManager::writePhoBranches(int info) {
    if (info == 1) {
        tree_->Branch("g_p",          &g_p_,           "g_p/D");
        tree_->Branch("g_theta",      &g_theta_,       "g_theta/D");
        tree_->Branch("g_phi",        &g_phi_,         "g_phi/D");
        tree_->Branch("g_sector",     &g_sector_,      "g_sector/I");
        tree_->Branch("g_vz",         &g_vz_,          "g_vz/D");
        tree_->Branch("g_status",     &g_status_,      "g_status/I");
        tree_->Branch("g_chi2pid",    &g_chi2pid_,     "g_chi2pid/D");
        tree_->Branch("detPho",       &detPho_,        "detPho/I");
        tree_->Branch("g_uPCAL",      &g_uPCAL_,       "g_uPCAL/D");
        tree_->Branch("g_vPCAL",      &g_vPCAL_,       "g_vPCAL/D");
        tree_->Branch("g_wPCAL",      &g_wPCAL_,       "g_wPCAL/D");
        tree_->Branch("g_uECIN",      &g_uECIN_,       "g_uECIN/D");
        tree_->Branch("g_vECIN",      &g_vECIN_,       "g_vECIN/D");
        tree_->Branch("g_wECIN",      &g_wECIN_,       "g_wECIN/D");
        tree_->Branch("g_uECOUT",     &g_uECOUT_,      "g_uECOUT/D");
        tree_->Branch("g_vECOUT",     &g_vECOUT_,      "g_vECOUT/D");
        tree_->Branch("g_wECOUT",     &g_wECOUT_,      "g_wECOUT/D");
        tree_->Branch("g_E_PCAL",     &g_E_PCAL_,      "g_E_PCAL/D");
        tree_->Branch("g_E_ECIN",     &g_E_ECIN_,      "g_E_ECIN/D");
        tree_->Branch("g_E_ECOUT",    &g_E_ECOUT_,     "g_E_ECOUT/D");
    }
}

// creates EPPI0 pi0 kinematics & exclusivity branches
void ProcessManager::writeEPPI0Branches() {
    // Pi0:
    tree_->Branch("pi0_p",        &pi0_p_,         "pi0_p/D");
    tree_->Branch("pi0_theta",    &pi0_theta_,     "pi0_theta/D");
    tree_->Branch("pi0_phi",      &pi0_phi_,       "pi0_phi/D");
    tree_->Branch("m_gg",         &m_gg_,          "m_gg/D");
    tree_->Branch("detPi0",       &detPi0_,        "detPi0/I");
    tree_->Branch("deltaPhi",     &deltaPhi_,      "deltaPhi/D");

    // Exclusivity branches:
    tree_->Branch("m2_miss",      &m2_miss_,       "m2_miss/D");
    tree_->Branch("m2_epX",       &m2_epX_,        "m2_epX/D");
    tree_->Branch("m2_epi0X",     &m2_epi0X_,      "m2_epi0X/D");
    tree_->Branch("E_miss",       &E_miss_,        "E_miss/D");
    tree_->Branch("E2_miss",      &E2_miss_,       "E2_miss/D");
    tree_->Branch("px_miss",      &px_miss_,       "px_miss/D");
    tree_->Branch("py_miss",      &py_miss_,       "py_miss/D");
    tree_->Branch("pz_miss",      &pz_miss_,       "pz_miss/D");
    tree_->Branch("pT_miss",      &pT_miss_,       "pT_miss/D");
    tree_->Branch("trentoPhi",    &trentoPhi_,     "trentoPhi/D");
}

// initializes ROOT Tree and links branches with corresponding event-time variables. NOTE: requires channel_ & topology_
void ProcessManager::rootTree() {
    tree_ = new TTree("Events", "CLAS12 event data");
    tree_->SetAutoSave(0);
    tree_->SetAutoFlush(5000);

    if (channel_ == "gen") {
        tree_->Branch("pid",       &pid_,           "pid/I");
        tree_->Branch("Q2_gen",    &Q2_gen_,        "Q2_gen/D");
        tree_->Branch("nu_gen",    &nu_gen_,        "nu_gen/D");
        tree_->Branch("Xb_gen",    &Xb_gen_,        "Xb_gen/D");
        tree_->Branch("y_gen",     &y_gen_,         "y_gen/D");
        tree_->Branch("W_gen",     &W_gen_,         "W_gen/D");
        tree_->Branch("t_gen",     &t_gen_,         "t_gen/D");
        tree_->Branch("e_pgen",    &e_pgen_,        "e_pgen/D");
        tree_->Branch("p_pgen",    &p_pgen_,        "p_pgen/D");
        writeEleBranches();
        return;
    }

    // General branches:
    tree_->Branch("runNum",        &runNum_,        "runNum/I");
    tree_->Branch("eventNum",      &eventNum_,      "eventNum/I");
    tree_->Branch("helicity",      &helicity_,      "helicity/I");
    tree_->Branch("pid",           &pid_,           "pid/I");
    tree_->Branch("charge",        &charge_,        "charge/I");
    tree_->Branch("chi2pid",       &chi2pid_,       "chi2pid/D");

    // Particle branches:
    writeEleBranches();
    writeProBranches();
    writePhoBranches();

    if (channel_ == "inclusiveGEMC") {
        tree_->Branch("Q2_gen",    &Q2_gen_,        "Q2_gen/D");
        tree_->Branch("nu_gen",    &nu_gen_,        "nu_gen/D");
        tree_->Branch("Xb_gen",    &Xb_gen_,        "Xb_gen/D");
        tree_->Branch("y_gen",     &y_gen_,         "y_gen/D");
        tree_->Branch("W_gen",     &W_gen_,         "W_gen/D");
        tree_->Branch("t_gen",     &t_gen_,         "t_gen/D");
        tree_->Branch("e_pgen",    &e_pgen_,        "e_pgen/D");
        tree_->Branch("p_pgen",    &p_pgen_,        "p_pgen/D");
    }

    else if (channel_ == "eppi0") {
        writeEPPI0Branches();
    }

}

// assigns vars depending on particle p, called by ProcessEvent on a particle basis
void ProcessManager::fillRecVars(clas12::region_particle* p, int ele_info, int pro_info, int pho_info) {

    int pid  = p->getPid();
    int detP = getDetector(p->par()->getStatus());

    pid_     = pid;
    charge_  = p->par()->getCharge();
    chi2pid_ = p->par()->getChi2Pid();

    if (ele_info > 0) {
        e_p_            = pid != 11 ? NAN   : p->getP();
        e_beta_         = pid != 11 ? NAN   : p->par()->getBeta();
        e_theta_        = pid != 11 ? NAN   : p->getTheta();
        e_phi_          = pid != 11 ? NAN   : p->getPhi();
        e_sector_       = pid != 11 ? -999  : p->getSector();
        e_vz_           = pid != 11 ? NAN   : p->par()->getVz();
        e_status_       = pid != 11 ? NAN   : p->par()->getStatus();
        detEle_         = pid != 11 ? -999  : detP;

        // DIS calculations:
        double e_Ef = std::sqrt(ELECTRON_MASS * ELECTRON_MASS + p->getP() * p->getP());
        TLorentzVector scatteredElectron(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(), e_Ef);
        Kinematics DIS(scatteredElectron, Ebeam_, PROTON_MASS);
        Q2_             = DIS.Q2();
        nu_             = DIS.nu();
        Xb_             = DIS.Xb();
        y_              = DIS.y();
        W_              = DIS.W();

        // Misc.
        e_xDC1_         = pid != 11 ? NAN : (detP == 1 && p->traj(DC, 6))  ? p->traj(DC, 6)->getX()  : NAN;
        e_yDC1_         = pid != 11 ? NAN : (detP == 1 && p->traj(DC, 6))  ? p->traj(DC, 6)->getY()  : NAN;
        e_xDC2_         = pid != 11 ? NAN : (detP == 1 && p->traj(DC, 6))  ? p->traj(DC, 6)->getX()  : NAN;
        e_yDC2_         = pid != 11 ? NAN : (detP == 1 && p->traj(DC, 18)) ? p->traj(DC, 18)->getY() : NAN;
        e_xDC3_         = pid != 11 ? NAN : (detP == 1 && p->traj(DC, 36)) ? p->traj(DC, 36)->getX() : NAN;
        e_yDC3_         = pid != 11 ? NAN : (detP == 1 && p->traj(DC, 36)) ? p->traj(DC, 36)->getY() : NAN;

        e_uPCAL_        = pid != 11 ? NAN : (detP == 1 && p->cal(1)) ? p->cal(1)->getLu()     : NAN;
        e_vPCAL_        = pid != 11 ? NAN : (detP == 1 && p->cal(1)) ? p->cal(1)->getLv()     : NAN;
        e_wPCAL_        = pid != 11 ? NAN : (detP == 1 && p->cal(1)) ? p->cal(1)->getLw()     : NAN;
        e_E_PCAL_       = pid != 11 ? NAN : (detP == 1 && p->cal(1)) ? p->cal(1)->getEnergy() : NAN;
    
        e_uECIN_        = pid != 11 ? NAN : (detP == 1 && p->cal(4)) ? p->cal(4)->getLu()     : NAN;
        e_vECIN_        = pid != 11 ? NAN : (detP == 1 && p->cal(4)) ? p->cal(4)->getLv()     : NAN;
        e_wECIN_        = pid != 11 ? NAN : (detP == 1 && p->cal(4)) ? p->cal(4)->getLw()     : NAN;
        e_E_ECIN_       = pid != 11 ? NAN : (detP == 1 && p->cal(4)) ? p->cal(4)->getEnergy() : NAN;

        e_uECOUT_       = pid != 11 ? NAN : (detP == 1 && p->cal(7)) ? p->cal(7)->getLu()     : NAN;
        e_vECOUT_       = pid != 11 ? NAN : (detP == 1 && p->cal(7)) ? p->cal(7)->getLv()     : NAN;
        e_wECOUT_       = pid != 11 ? NAN : (detP == 1 && p->cal(7)) ? p->cal(7)->getLw()     : NAN;
        e_E_ECOUT_      = pid != 11 ? NAN : (detP == 1 && p->cal(7)) ? p->cal(7)->getEnergy() : NAN;

    }

    if (pro_info > 0) {

        p_p_            = pid != 2212 ? NAN   : p->getP();
        p_beta_         = pid != 2212 ? NAN   : p->par()->getBeta();
        p_theta_        = pid != 2212 ? NAN   : p->getTheta();
        p_phi_          = pid != 2212 ? NAN   : p->getPhi();
        p_sector_       = pid != 2212 ? -999  : p->getSector();
        p_vz_           = pid != 2212 ? NAN   : p->par()->getVz();
        p_status_       = pid != 2212 ? NAN   : p->par()->getStatus();
        detPro_         = pid != 2212 ? -999  : detP;

        double p_Ef = std::sqrt(PROTON_MASS * PROTON_MASS + p->getP() * p->getP());
        TLorentzVector target(0,0,0,PROTON_MASS);
        TLorentzVector finalProton(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(), p_Ef);
        TLorentzVector t = finalProton - target;
        t_              = pid != 2212 ? NAN : -1 * t.Mag2();

        p_xDC1_         = pid != 2212 ? NAN : (detP == 1 && p->traj(DC, 6))   ? p->traj(DC, 6)->getX()  : NAN;
        p_yDC1_         = pid != 2212 ? NAN : (detP == 1 && p->traj(DC, 6))   ? p->traj(DC, 6)->getY()  : NAN;
        p_xDC2_         = pid != 2212 ? NAN : (detP == 1 && p->traj(DC, 18))  ? p->traj(DC, 18)->getX() : NAN;
        p_yDC2_         = pid != 2212 ? NAN : (detP == 1 && p->traj(DC, 18))  ? p->traj(DC, 18)->getY() : NAN;
        p_xDC3_         = pid != 2212 ? NAN : (detP == 1 && p->traj(DC, 36))  ? p->traj(DC, 36)->getX() : NAN;
        p_yDC3_         = pid != 2212 ? NAN : (detP == 1 && p->traj(DC, 36))  ? p->traj(DC, 36)->getY() : NAN;

        p_edge_cvt1_    = pid != 2212 ? NAN : (detP == 2 && p->traj(CVT, 1))  ? p->traj(CVT, 1)->getEdge()  : NAN;
        p_edge_cvt3_    = pid != 2212 ? NAN : (detP == 2 && p->traj(CVT, 3))  ? p->traj(CVT, 3)->getEdge()  : NAN;
        p_edge_cvt5_    = pid != 2212 ? NAN : (detP == 2 && p->traj(CVT, 5))  ? p->traj(CVT, 5)->getEdge()  : NAN;
        p_edge_cvt7_    = pid != 2212 ? NAN : (detP == 2 && p->traj(CVT, 7))  ? p->traj(CVT, 7)->getEdge()  : NAN;
        p_edge_cvt12_   = pid != 2212 ? NAN : (detP == 2 && p->traj(CVT, 12)) ? p->traj(CVT, 12)->getEdge() : NAN;

        double x_cvt    = pid != 2212 ? NAN : (detP == 2 && p->traj(CVT, 1))  ? p->traj(CVT, 1)->getX() : NAN;
        double y_cvt    = pid != 2212 ? NAN : (detP == 2 && p->traj(CVT, 1))  ? p->traj(CVT, 1)->getY() : NAN;
        double z_cvt    = pid != 2212 ? NAN : (detP == 2 && p->traj(CVT, 1))  ? p->traj(CVT, 1)->getZ() : NAN;
        double r_cvt    = pid != 2212 ? NAN : std::sqrt(x_cvt * x_cvt + y_cvt * y_cvt + z_cvt * z_cvt);

        p_theta_cvt_    = pid != 2212 ? NAN : (detP == 2 && !isnan(r_cvt)) ? std::acos(z_cvt / r_cvt) : NAN;
        p_phi_cvt_      = pid != 2212 ? NAN : (detP == 2 && !isnan(x_cvt) && !isnan(y_cvt)) ? std::atan2(y_cvt, x_cvt) : NAN;

    }

    if (pho_info > 0) {

        g_p_            = pid != 22 ? NAN   : p->getP();
        g_theta_        = pid != 22 ? NAN   : p->getTheta();
        g_phi_          = pid != 22 ? NAN   : p->getPhi();
        g_sector_       = pid != 22 ? -999  : p->getSector();
        g_vz_           = pid != 22 ? NAN   : p->par()->getVz();
        g_status_       = pid != 22 ? NAN   : p->par()->getStatus();
        detPho_         = pid != 22 ? -999  : detP;

        g_uPCAL_        = pid != 22 ? NAN : (detP == 1 && p->cal(1)) ? p->cal(1)->getLu()     : NAN;
        g_vPCAL_        = pid != 22 ? NAN : (detP == 1 && p->cal(1)) ? p->cal(1)->getLv()     : NAN;
        g_wPCAL_        = pid != 22 ? NAN : (detP == 1 && p->cal(1)) ? p->cal(1)->getLw()     : NAN;
        g_E_PCAL_       = pid != 22 ? NAN : (detP == 1 && p->cal(1)) ? p->cal(1)->getEnergy() : NAN;
    
        g_uECIN_        = pid != 22 ? NAN : (detP == 1 && p->cal(4)) ? p->cal(4)->getLu()     : NAN;
        g_vECIN_        = pid != 22 ? NAN : (detP == 1 && p->cal(4)) ? p->cal(4)->getLv()     : NAN;
        g_wECIN_        = pid != 22 ? NAN : (detP == 1 && p->cal(4)) ? p->cal(4)->getLw()     : NAN;
        g_E_ECIN_       = pid != 22 ? NAN : (detP == 1 && p->cal(4)) ? p->cal(4)->getEnergy() : NAN;

        g_uECOUT_       = pid != 22 ? NAN : (detP == 1 && p->cal(7)) ? p->cal(7)->getLu()     : NAN;
        g_vECOUT_       = pid != 22 ? NAN : (detP == 1 && p->cal(7)) ? p->cal(7)->getLv()     : NAN;
        g_wECOUT_       = pid != 22 ? NAN : (detP == 1 && p->cal(7)) ? p->cal(7)->getLw()     : NAN;
        g_E_ECOUT_      = pid != 22 ? NAN : (detP == 1 && p->cal(7)) ? p->cal(7)->getEnergy() : NAN;

    }

}

// assigns electron vars from electron `ele`, called by ProcessEvent() on an event basis
void ProcessManager::fillEleVars(clas12::region_particle* ele, int info) {
    if (info > 0) {
        e_p_            = ele->getP();
        e_theta_        = ele->getTheta();
        e_beta_         = ele->par()->getBeta();
        e_phi_          = ele->getPhi();
        e_sector_       = ele->getSector();
        e_vz_           = ele->par()->getVz();
        e_status_       = ele->par()->getStatus();
        e_chi2pid_      = ele->par()->getChi2Pid();
        detEle_         = getDetector(ele->par()->getStatus());

        //DIS:
        double e_Ef = std::sqrt(ELECTRON_MASS * ELECTRON_MASS + ele->getP() * ele->getP());
        TLorentzVector scatteredElectron(ele->par()->getPx(), ele->par()->getPy(), ele->par()->getPz(), e_Ef);
        // std::cout << "PM initializes e_f with px = " << scatteredElectron.Px() << ", py = " << scatteredElectron.Py() <<
        //             ", pz = " << scatteredElectron.Pz() << ", E = " << scatteredElectron.E() << std::endl;
        Kinematics DIS(scatteredElectron, Ebeam_, PROTON_MASS);
        Q2_             = DIS.Q2();
        nu_             = DIS.nu();
        Xb_             = DIS.Xb();
        y_              = DIS.y();
        W_              = DIS.W();
    
        //std::cout << "PM about to fill Q2 = " << Q2_ << ", Xb_ = " << Xb_ << std::endl;

        e_xDC1_         = ele->traj(DC, 6)  ? ele->traj(DC, 6)->getX()  : NAN;
        e_yDC1_         = ele->traj(DC, 6)  ? ele->traj(DC, 6)->getY()  : NAN;
        e_xDC2_         = ele->traj(DC, 6)  ? ele->traj(DC, 6)->getX()  : NAN;
        e_yDC2_         = ele->traj(DC, 18) ? ele->traj(DC, 18)->getY() : NAN;
        e_xDC3_         = ele->traj(DC, 36) ? ele->traj(DC, 36)->getX() : NAN;
        e_yDC3_         = ele->traj(DC, 36) ? ele->traj(DC, 36)->getY() : NAN;

        e_uPCAL_        = ele->cal(1) ? ele->cal(1)->getLu()     : NAN;
        e_vPCAL_        = ele->cal(1) ? ele->cal(1)->getLv()     : NAN;
        e_wPCAL_        = ele->cal(1) ? ele->cal(1)->getLw()     : NAN;
        e_E_PCAL_       = ele->cal(1) ? ele->cal(1)->getEnergy() : NAN;
    
        e_uECIN_        = ele->cal(4) ? ele->cal(4)->getLu()     : NAN;
        e_vECIN_        = ele->cal(4) ? ele->cal(4)->getLv()     : NAN;
        e_wECIN_        = ele->cal(4) ? ele->cal(4)->getLw()     : NAN;
        e_E_ECIN_       = ele->cal(4) ? ele->cal(4)->getEnergy() : NAN;

        e_uECOUT_       = ele->cal(7) ? ele->cal(7)->getLu()     : NAN;
        e_vECOUT_       = ele->cal(7) ? ele->cal(7)->getLv()     : NAN;
        e_wECOUT_       = ele->cal(7) ? ele->cal(7)->getLw()     : NAN;
        e_E_ECOUT_      = ele->cal(7) ? ele->cal(7)->getEnergy() : NAN;
    }

}

// assigns proton vars from proton `pro`, called by ProcessEvent() on an event basis
void ProcessManager::fillProVars(clas12::region_particle* pro, int info) {
    if (info > 0) {

        p_p_            = pro->getP();
        p_beta_         = pro->par()->getBeta();
        p_theta_        = pro->getTheta();
        p_phi_          = pro->getPhi();
        p_sector_       = pro->getSector();
        p_vz_           = pro->par()->getVz();
        p_status_       = pro->par()->getStatus();
        p_chi2pid_      = pro->par()->getChi2Pid();
        detPro_         = getDetector(pro->par()->getStatus());

        double p_Ef = std::sqrt(PROTON_MASS * PROTON_MASS + pro->getP() * pro->getP());
        TLorentzVector target(0,0,0,PROTON_MASS);
        TLorentzVector finalProton(pro->par()->getPx(), pro->par()->getPy(), pro->par()->getPz(), p_Ef);
        TLorentzVector t = finalProton - target;
        t_              = -1 * t.Mag2();

        p_xDC1_         = (detPro_ == 1 && pro->traj(DC, 6))   ? pro->traj(DC, 6)->getX()  : NAN;
        p_yDC1_         = (detPro_ == 1 && pro->traj(DC, 6))   ? pro->traj(DC, 6)->getY()  : NAN;
        p_xDC2_         = (detPro_ == 1 && pro->traj(DC, 18))  ? pro->traj(DC, 18)->getX() : NAN;
        p_yDC2_         = (detPro_ == 1 && pro->traj(DC, 18))  ? pro->traj(DC, 18)->getY() : NAN;
        p_xDC3_         = (detPro_ == 1 && pro->traj(DC, 36))  ? pro->traj(DC, 36)->getX() : NAN;
        p_yDC3_         = (detPro_ == 1 && pro->traj(DC, 36))  ? pro->traj(DC, 36)->getY() : NAN;

        p_edge_cvt1_    = (detPro_ == 2 && pro->traj(CVT, 1))  ? pro->traj(CVT, 1)->getEdge()  : NAN;
        p_edge_cvt3_    = (detPro_ == 2 && pro->traj(CVT, 3))  ? pro->traj(CVT, 3)->getEdge()  : NAN;
        p_edge_cvt5_    = (detPro_ == 2 && pro->traj(CVT, 5))  ? pro->traj(CVT, 5)->getEdge()  : NAN;
        p_edge_cvt7_    = (detPro_ == 2 && pro->traj(CVT, 7))  ? pro->traj(CVT, 7)->getEdge()  : NAN;
        p_edge_cvt12_   = (detPro_ == 2 && pro->traj(CVT, 12)) ? pro->traj(CVT, 12)->getEdge() : NAN;

        double x_cvt    = (detPro_ == 2 && pro->traj(CVT, 1))  ? pro->traj(CVT, 1)->getX() : NAN;
        double y_cvt    = (detPro_ == 2 && pro->traj(CVT, 1))  ? pro->traj(CVT, 1)->getY() : NAN;
        double z_cvt    = (detPro_ == 2 && pro->traj(CVT, 1))  ? pro->traj(CVT, 1)->getZ() : NAN;
        double r_cvt    = std::sqrt(x_cvt * x_cvt + y_cvt * y_cvt + z_cvt * z_cvt);

        p_theta_cvt_    =  (detPro_ == 2 && !isnan(r_cvt)) ? std::acos(z_cvt / r_cvt) : NAN;
        p_phi_cvt_      =  (detPro_ == 2 && !isnan(x_cvt) && !isnan(y_cvt)) ? std::atan2(y_cvt, x_cvt) : NAN;
    }
}

// assigns photon vars from photon `pho`, called by ProcessEvent() on an event basis
void ProcessManager::fillPhoVars(clas12::region_particle* pho, int info) {
    if (info > 0) {
        g_p_            = pho->getP();
        g_theta_        = pho->getTheta();
        g_phi_          = pho->getPhi();
        g_sector_       = pho->getSector();
        g_vz_           = pho->par()->getVz();
        g_status_       = pho->par()->getStatus();
        g_chi2pid_      = pho->par()->getChi2Pid();
        detPho_         = getDetector(pho->par()->getStatus());

        g_uPCAL_        = pho->cal(1) ? pho->cal(1)->getLu()     : NAN;
        g_vPCAL_        = pho->cal(1) ? pho->cal(1)->getLv()     : NAN;
        g_wPCAL_        = pho->cal(1) ? pho->cal(1)->getLw()     : NAN;
        g_E_PCAL_       = pho->cal(1) ? pho->cal(1)->getEnergy() : NAN;
    
        g_uECIN_        = pho->cal(4) ? pho->cal(4)->getLu()     : NAN;
        g_vECIN_        = pho->cal(4) ? pho->cal(4)->getLv()     : NAN;
        g_wECIN_        = pho->cal(4) ? pho->cal(4)->getLw()     : NAN;
        g_E_ECIN_       = pho->cal(4) ? pho->cal(4)->getEnergy() : NAN;

        g_uECOUT_       = pho->cal(7) ? pho->cal(7)->getLu()     : NAN;
        g_vECOUT_       = pho->cal(7) ? pho->cal(7)->getLv()     : NAN;
        g_wECOUT_       = pho->cal(7) ? pho->cal(7)->getLw()     : NAN;
        g_E_ECOUT_      = pho->cal(7) ? pho->cal(7)->getEnergy() : NAN;
    }
}

// assigns EPPI0 exclusivity vars from Kinematics instance `Kin`
void ProcessManager::fillEPPI0Vars(Kinematics Kin) {
        m2_miss_   = Kin.missingP4().M2();
        m2_epX_    = Kin.lv_epX().M2();
        m2_epi0X_  = Kin.lv_epi0X().M2();
        E_miss_    = Kin.missingP4().E();
        E2_miss_   = Kin.missingP4().E() * Kin.missingP4().E();
        px_miss_   = Kin.missingP4().Px();
        py_miss_   = Kin.missingP4().Py();
        pz_miss_   = Kin.missingP4().Pz();
        pT_miss_   = Kin.missingP4().Pt();
        trentoPhi_ = Kin.phiT();
}

// writes ROOT Tree to `outFile_` and frees up memory allocated to `tree_` and `outFile`. Called AFTER processing
void ProcessManager::finalize(const std::string& output_file) {
    outFile_ = new TFile(output_file.c_str(), "RECREATE");
    if (outFile_) {
        outFile_->cd();
        if (tree_) {
            tree_->Write();
            delete tree_;
            tree_ = nullptr;
        }
        outFile_->Close();
        delete outFile_;
        outFile_ = nullptr;
    } else if (tree_) {
        // Tree exists but no outfile? Just delete tree to avoid leaks.
        delete tree_;
        tree_ = nullptr;
    }
}

// MAIN WORKFLOW PERFORMED HERE:
void ProcessManager::processEvent(clas12::clas12reader& c12) {
    if (channel_ == "trigger") {
        auto electrons = c12.getByID(11);
        if (electrons.size() < 1) return;
        auto& trigEle = electrons[0];
        int detP = getDetector(trigEle->par()->getStatus());
        bool keepEle = ((detP == 0 && fiducialCuts_->passesFT(trigEle)) || 
                                (detP == 1 && fiducialCuts_->passesDC(trigEle, torus_) && fiducialCuts_->passesECAL(trigEle)));
        if (!keepEle) return;

        runNum_   = c12.runconfig()->getRun();
        eventNum_ = c12.runconfig()->getEvent();
        helicity_ = c12.event()->getHelicity();

        fillEleVars(trigEle);
        tree_->Fill();
        numFills_++;
    }
    else if (channel_ == "gen") {

        // Generated particles
        auto mcParticles = c12.mcparts();
        int numGen = mcParticles->getRows();       

        for (int j = 0; j < numGen; j++) {
            int gen_pid = mcParticles->getPid(j);
            pid_ = gen_pid;
    
            TVector3 vec_gen_p(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j));

            e_pgen_ = p_pgen_ = NAN;
            Q2_gen_ = nu_gen_ = Xb_gen_ = y_gen_ = W_gen_ = t_gen_ = NAN;

            if (gen_pid == 11) {

                e_pgen_   = vec_gen_p.Mag();
                
                // DIS calculations:
                double e_Ef_gen = std::sqrt(ELECTRON_MASS * ELECTRON_MASS + vec_gen_p.Mag() * vec_gen_p.Mag());
                TLorentzVector scatteredGenElectron(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j), e_Ef_gen);
                Kinematics GEN_DIS(scatteredGenElectron, Ebeam_, PROTON_MASS);
                Q2_gen_  = GEN_DIS.Q2();
                nu_gen_  = GEN_DIS.nu();
                Xb_gen_  = GEN_DIS.Xb();
                y_gen_   = GEN_DIS.y();
                W_gen_   = GEN_DIS.W();
            } 

            else if (gen_pid == 2212) {

                p_pgen_   = vec_gen_p.Mag();

                double p_Ef_gen = std::sqrt(PROTON_MASS * PROTON_MASS + vec_gen_p.Mag() * vec_gen_p.Mag());
                TLorentzVector target(0,0,0,PROTON_MASS);
                TLorentzVector finalProton(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j), p_Ef_gen);
                TLorentzVector t = finalProton - target;
                t_gen_ = -1 * t.Mag2();
            }

            tree_->Fill();
            numFills_++;
        }

    }
    // inclusive channel: accepts electrons, protons, photons
    else if (channel_ == "inclusive") {
        
        auto recParticles = c12.getDetParticles();   
        int numRec = recParticles.size();            

        for (int i = 0; i < numRec; i++) {
            auto& p = recParticles[i];

            int pid = p->getPid();
            int detP = getDetector(p->par()->getStatus());

            if (pid == 11) {
                bool keepEle = ((detP == 0 && fiducialCuts_->passesFT(p)) || 
                                (detP == 1 && fiducialCuts_->passesDC(p, torus_) && fiducialCuts_->passesECAL(p)));
                if (!keepEle) return;
            } 

            else if (pid == 2212) {
                bool keepPro = ((detP == 1 && fiducialCuts_->passesDC(p, torus_)) || 
                                (detP == 2 && fiducialCuts_->passesCVT(p)));
                if (!keepPro) return;
            }

            else if (pid == 22) {
                bool keepPho = (detP == 0 && fiducialCuts_->passesFT(p)) || (detP == 1 && fiducialCuts_->passesECAL(p));
                if (!keepPho) return;
            }

            runNum_   = c12.runconfig()->getRun();
            eventNum_ = c12.runconfig()->getEvent();
            helicity_ = c12.event()->getHelicity();

            fillRecVars(p);
            tree_->Fill();
            numFills_++;
        }
    }
    else if (channel_ == "inclusiveGEMC") {
        
        // Reconstructed particles
        auto recParticles = c12.getDetParticles();
        // Generated particles
        auto mcParticles = c12.mcparts();

        int numGen = mcParticles->getRows();  
        int numRec = recParticles.size();     

        for (int j = 0; j < numGen; j++) {
            int gen_pid = mcParticles->getPid(j);
            TVector3 vec_gen_p(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j));
            double bestDeltaTheta = 999.0;
            int bestMatchIndex = -1;

            for (int i = 0; i < numRec; i++) {
                auto& rec = recParticles[i];
                if (rec->par()->getPid() != gen_pid) continue;

                TVector3 vec_rec_p(rec->par()->getPx(), rec->par()->getPy(), rec->par()->getPz());
                // double deltaPhi = fabs(vec_gen_p.Phi() - vec_rec_p.Phi()) * 180/M_PI;
                // double deltaTheta = fabs(vec_gen_p.Theta() - vec_rec_p.Theta()) * 180/M_PI;
                double deltaTheta = vec_gen_p.Angle(vec_rec_p) * 180/M_PI; // in degrees

                if (deltaTheta < bestDeltaTheta) {
                    bestDeltaTheta = deltaTheta;
                    bestMatchIndex = i;
                }
            }

            // Threshold for match quality:
            if (bestMatchIndex < 0 || bestDeltaTheta > 3) continue;

            eventNum_ = c12.runconfig()->getEvent();
            auto& p   = recParticles[bestMatchIndex];
            int detP  = getDetector(p->par()->getStatus());

            e_pgen_ = p_pgen_ = NAN;
            Q2_gen_ = nu_gen_ = Xb_gen_ = y_gen_ = W_gen_ = t_gen_ = NAN;

            if (gen_pid == 11) {
                bool keepEle = passesVertexCut(p) && ((detP == 0 && fiducialCuts_->passesFT(p)) || 
                                (detP == 1 && fiducialCuts_->passesDC(p, torus_) && fiducialCuts_->passesECAL(p)));
                if (!keepEle) return;

                e_pgen_   = vec_gen_p.Mag();
                
                // DIS calculations:
                double e_Ef_gen = std::sqrt(ELECTRON_MASS * ELECTRON_MASS + vec_gen_p.Mag() * vec_gen_p.Mag());
                TLorentzVector scatteredGenElectron(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j), e_Ef_gen);
                Kinematics GEN_DIS(scatteredGenElectron, Ebeam_, PROTON_MASS);
                Q2_gen_  = GEN_DIS.Q2();
                nu_gen_  = GEN_DIS.nu();
                Xb_gen_  = GEN_DIS.Xb();
                y_gen_   = GEN_DIS.y();
                W_gen_   = GEN_DIS.W();
            } 

            else if (gen_pid == 2212) {
                bool keepPro = passesVertexCut(p) && ((detP == 1 && fiducialCuts_->passesDC(p, torus_)) || 
                                (detP == 2 && fiducialCuts_->passesCVT(p)));
                if (!keepPro) return;

                p_pgen_   = vec_gen_p.Mag();

                double p_Ef_gen = std::sqrt(PROTON_MASS * PROTON_MASS + vec_gen_p.Mag() * vec_gen_p.Mag());
                TLorentzVector target(0,0,0,PROTON_MASS);
                TLorentzVector finalProton(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j), p_Ef_gen);
                TLorentzVector t = finalProton - target;
                t_gen_ = -1 * t.Mag2();
            }

            else if (gen_pid == 22) {
                bool keepPho = (detP == 0 && fiducialCuts_->passesFT(p)) || (detP == 1 && fiducialCuts_->passesECAL(p));
                if (!keepPho) return;
            }

            else return;

            fillRecVars(p);
            tree_->Fill();
            numFills_++;
        }
    }
    else if (channel_ == "matchingWIP") {

        auto recParticles = c12.getDetParticles();
        auto mcParticles = c12.mcparts();

        int numGen = mcParticles->getRows();  // number of MC particles
        int numRec = recParticles.size();     // number of reconstructed particles

        for (int i = 0; i < numGen; i++) {
            mcParticles->setEntry(i);  

            if (!mcParticles->isMatched()) {
                std::cout << "MC particle " << i << " was NOT reconstructed.\n";
                continue;
            }

            int recIndex = mcParticles->getMatch()->getPindex();  // index of matched reconstructed particle

            auto match = mcParticles->getMatch();
            if (!match || !mcParticles->isMatched()) {
                std::cout << "MC particle " << i << " has no valid match.\n";
                continue;


                std::cout << "MC PID: " << mcParticles->getPid()
                        << " matched to REC PID: " << recParticles[recIndex]->par()->getPid()
                        << " with pGEN: " << mcParticles->getP()
                        << ", pREC: " << recParticles[recIndex]->par()->getP()
                        << std::endl;
            }

        } 

        for(auto p : c12.getDetParticles()){
            if(p->mc()->isMatched()){ //this particle has an mc match
                if(p->mc()->getMatch()->getQuality()>0.98){

                    std::cout<<p->par()->getEntry()<<" rec pid "<<p->par()->getPid()<<", mc pid "<<p->mc()->getPid()<<", mcindex "<< p->mc()->getMatch()->getMCindex()<<", quality = " << p->mc()->getMatch()->getQuality() << std::endl;
                }
            }
        }

    }
    else if (channel_ == "eppi0") {

        auto electrons = c12.getByID(11);   // PID 11 = electron
        auto protons   = c12.getByID(2212); // PID 2212 = proton
        auto photons   = c12.getByID(22);   // PID 22 = photon

        // Reject if not EPPI0 final state, according to eventbuilder:
        if (photons.size() != 2 || electrons.size() != 1 || protons.size() != 1) return;

        clas12::region_particle* best_electron = nullptr;
        clas12::region_particle* best_proton   = nullptr;
        clas12::region_particle* best_g1       = nullptr;
        clas12::region_particle* best_g2       = nullptr;
        
        TLorentzVector bestPi0;
        int detPi0;

        double  leastDeltaM = 999;
        double  m_bestPi0 = 999;

        for (const auto& ele : electrons) {

            // Reject if electron vertex isn't in target fiducial volume
            if (!passesVertexCut(ele)) continue;

            int detEle = getDetector(ele->par()->getStatus());
            // Reject if electron fails fiducial cuts (ONLY if cuts are listed!)
            if (detEle == 0 && !fiducialCuts_->passesFT(ele)) continue;
            if (detEle == 1 && (!fiducialCuts_->passesDC(ele, torus_) || !fiducialCuts_->passesECAL(ele))) continue;

            for (const auto& prot : protons) {

                // Reject if proton vertex isn't in target fiducial volume
                if (!passesVertexCut(prot)) continue;

                // Reject if topology has been specified (underscored private variables) && proton fails the topology
                int detPro = getDetector(prot->par()->getStatus());
                if (requireTopology_ && detPro != detPro_) continue;

                // Reject if proton fails fiducial cuts (ONLY if cuts are listed!)
                if (detPro == 1 && !fiducialCuts_->passesDC(prot, torus_)) continue;
                if (detPro == 2 && !fiducialCuts_->passesCVT(prot)) continue;

                for (size_t i = 0; i < photons.size(); ++i) {
                    for (size_t j = i + 1; j < photons.size(); ++j) {

                        auto& g1 = photons[i];
                        auto& g2 = photons[j];

                        int det1 = getDetector(g1->par()->getStatus());
                        int det2 = getDetector(g2->par()->getStatus());

                        // Reject if either photon is detected in CD
                        if (det1 == 2 || det2 == 2) continue;

                        int sec1 = g1->getSector();
                        int sec2 = g2->getSector();

                        // Reject if photons are not detected in same sector
                        if (sec1 != sec2) continue;

                        // Reject if topology has been specified && either photon fails the topology
                        if (requireTopology_ && (det1 != detPho_ || det2 != detPho_)) continue;

                        // Reject if either photon is detected in FT and fails FT fiducial cuts (fails ONLY if FTstandardCut is listed!)
                        if ((det1 == 0 && !fiducialCuts_->passesFT(g1)) || (det2 == 0 && !fiducialCuts_->passesFT(g2))) continue;

                        //Reject if either photon is detected in FD and fails ECAL fiducial cuts (fails ONLY if ECAL cuts are listed!)
                        if ((det1 == 1 && !fiducialCuts_->passesECAL(g1)) || (det2 == 1 && !fiducialCuts_->passesECAL(g2))) continue;

                        TLorentzVector lv_g1, lv_g2;
                        lv_g1.SetXYZM(g1->par()->getPx(), g1->par()->getPy(), g1->par()->getPz(), 0.0);
                        lv_g2.SetXYZM(g2->par()->getPx(), g2->par()->getPy(), g2->par()->getPz(), 0.0);
                        TLorentzVector candidatePi0 = lv_g1 + lv_g2;

                        double deltaM = std::abs(candidatePi0.M() - PI0_MASS);
                        if (deltaM > 0.025) continue;

                        if (deltaM < leastDeltaM) {
                            m_bestPi0      = candidatePi0.M();
                            detPi0         = det1; // Assign to pi0 the detector of the photon
                            leastDeltaM    = deltaM;
                            bestPi0        = candidatePi0;
                            best_electron  = ele;
                            best_proton    = prot;
                            best_g1        = g1;
                            best_g2        = g2;
                        }
                    }
                }
            }
        }

        // Reject if no suitable EPPI0 candidate has been found, i.e., best_electron == nullptr:
        if (!best_electron) return; 

        runNum_         = c12.runconfig()->getRun();
        eventNum_       = c12.runconfig()->getEvent();
        helicity_       = c12.event()->getHelicity();

        // ELECTRON INFO:
        fillEleVars(best_electron);

        // PROTON INFO:
        fillProVars(best_proton);

        // PHOTON INFO:
        fillPhoVars(best_g1);

        // PION INFO:
        pi0_p_          = bestPi0.P();
        pi0_theta_      = bestPi0.Theta();
        pi0_phi_        = bestPi0.Phi();
        m_gg_           = m_bestPi0;
        detPi0_         = detPi0;

        // DIS + EPPI0 INFO:
        double e_Ef = std::sqrt(ELECTRON_MASS * ELECTRON_MASS + best_electron->getP() * best_electron->getP());
        double p_Ef = std::sqrt(PROTON_MASS   * PROTON_MASS   + best_proton->getP()   * best_proton->getP());

        TLorentzVector scatteredElectron(best_electron->par()->getPx(), best_electron->par()->getPy(), 
            best_electron->par()->getPz(), e_Ef);

        TLorentzVector finalProton(best_proton->par()->getPx(), best_proton->par()->getPy(),
            best_proton->par()->getPz(), p_Ef);

        Kinematics EPPI0(scatteredElectron, finalProton, bestPi0, Ebeam_);

        // Reject if Q2 < 1, W < 2, or y > 0.8, i.e., not in standard DIS region
        if (!EPPI0.channelCheck()) return;

        fillEPPI0Vars(EPPI0);
    
        deltaPhi_  = bestPi0.Phi() - EPPI0.lv_epX().Phi();

        // FILL TREE, ONCE PER EVENT:
        tree_->Fill();
        numFills_++;
    }
    eventsProcessed_++;
}

