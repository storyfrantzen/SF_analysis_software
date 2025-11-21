#pragma once

#include "TObject.h"
#include "TLorentzVector.h"
#include "clas12reader.h"
#include "PhysicalConstants.h"
#include <string>
#include <unordered_set>
#include <cmath> 

using namespace clas12;

struct EventVars : public TObject {
    // ALWAYS FILLED: //
    int runNum, eventNum, helicity;

    double charge;

    // FILLED IF REQUESTED: //
    int numElectrons, numProtons, numPhotons;

    void flush() {
        runNum = eventNum = helicity = numElectrons = numProtons = numPhotons = -9999;
        charge = NAN;
    }
    void fill(const std::unordered_set<std::string>& enabledBranches, clas12::clas12reader& c12);
    ClassDef(EventVars, 1);
};
struct RecVars : public TObject {

    // ALWAYS FILLED: //
    int    pid, charge, status, det, sector;
    double p, beta, theta, phi, chi2pid;

    // FILLED IF REQUESTED: //
    double px, py, pz, vx, vy, vz, time; 
    double xFT, yFT, xDC1, yDC1, xDC2, yDC2, xDC3, yDC3;
    double xPCAL, yPCAL, uPCAL, vPCAL, wPCAL, uECIN, vECIN, wECIN, uECOUT, vECOUT, wECOUT;
    double E_PCAL, E_ECIN, E_ECOUT;
    double edge_cvt1, edge_cvt3, edge_cvt5, edge_cvt7, edge_cvt12, theta_cvt, phi_cvt;

    void flush() {
        pid = charge = status = det = sector = -9999;
        p = beta = theta = phi = chi2pid = NAN;
        px = py = pz = vx = vy = vz = time = NAN;
        xFT = yFT = xDC1 = yDC1 = xDC2 = yDC2 = xDC3 = yDC3 = NAN;
        xPCAL = yPCAL = uPCAL = vPCAL = wPCAL = uECIN = vECIN = wECIN = uECOUT = vECOUT = wECOUT = NAN;
        E_PCAL = E_ECIN = E_ECOUT = NAN;
        edge_cvt1 = edge_cvt3 = edge_cvt5 = edge_cvt7 = edge_cvt12 = theta_cvt = phi_cvt = NAN;
    }
    void fill(const std::unordered_set<std::string>& enabledBranches, clas12::region_particle* part);
    void fill(const std::unordered_set<std::string>& enabledBranches, clas12::region_particle* rec,
                   double pIn, double thetaIn, double phiIn);
    ClassDef(RecVars, 1);
};
struct GenVars : public TObject {
    // ALWAYS FILLED: //
    int pid;
    double p, theta, phi;

    void flush() {
        pid = -9999;
        p = theta = phi = NAN;
    }
    void fill(clas12::mcparticle* mcParticles, int j);
    ClassDef(GenVars, 1);
};
struct DISVars : public TObject {
    // ALWAYS FILLED from trigger electron, redundant for successive event particles //
    double Q2, nu, Xb, y, W;

    void fill(const TLorentzVector& lv_ePrime, double ebeam, double m_target=PROTON_MASS);
    void flush() {
        Q2 = nu = Xb = y = W = NAN;
    };
    ClassDef(DISVars, 1);
};
struct EPPI0Vars : public TObject {
    // ALWAYS FILLED: //
    double pi0_p, pi0_theta, pi0_phi, pi0_deltaPhi, pi0_thetaX;
    double m_gg, m2_miss, m2_epX, m2_epi0X, m_eggX;
    double E_miss, px_miss, py_miss, pz_miss, pT_miss;
    double theta_e_g1, theta_e_g2, theta_g1_g2;
    double t, trentoPhi;
    
    void flush() {
        pi0_p = pi0_theta = pi0_phi = pi0_deltaPhi = NAN;
        m_gg = m2_miss = m2_epX = m2_epi0X = NAN;
        E_miss = px_miss = py_miss = pz_miss = pT_miss = NAN;
        theta_e_g1 = theta_e_g2 = theta_g1_g2 = NAN;
        t = trentoPhi = NAN;
    }

    void fill(const TLorentzVector& lv_ePrime, const TLorentzVector& lv_pPrime, 
              const TLorentzVector& lv_g1, const TLorentzVector& lv_g2, double ebeam, double m_target=PROTON_MASS);


    private:
    // Helper: compute Trento phi
    // lv_beam, lv_target, lv_ePrime: lab-frame vectors
    // lv_hadron: the hadron used to define the hadron plane (e.g., recoil proton or π0)
    double computeTrentoPhi(const TLorentzVector& lv_beam,
                            const TLorentzVector& lv_target,
                            const TLorentzVector& lv_ePrime,
                            const TLorentzVector& lv_hadron) const;
    ClassDef(EPPI0Vars, 1);

};