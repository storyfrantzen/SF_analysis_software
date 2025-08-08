#pragma once

#include "TObject.h"
#include "clas12reader.h"
#include <cmath> 

using namespace clas12;

struct EventVars : public TObject {
    int runNum, eventNum, helicity;

    void fill(clas12::clas12reader& c12);
    ClassDef(EventVars, 1);
};
struct RecVars : public TObject {
    int    pid, charge, status, det, sector;
    double p, beta, theta, phi, px, py, pz, vx, vy, vz, chi2pid, time; 
    double xFT, yFT, xDC1, yDC1, xDC2, yDC2, xDC3, yDC3;
    double xPCAL, yPCAL, uPCAL, vPCAL, wPCAL, uECIN, vECIN, wECIN, uECOUT, vECOUT, wECOUT;
    double E_PCAL, E_ECIN, E_ECOUT;
    double edge_cvt1, edge_cvt3, edge_cvt5, edge_cvt7, edge_cvt12, theta_cvt, phi_cvt;

    void fill(clas12::region_particle* part);
    ClassDef(RecVars, 1);
};
struct GenVars : public TObject {
    int pid;
    double p, theta, phi, Q2, nu, Xb, y, W, t;

    ClassDef(GenVars, 1);
};
struct DISVars : public TObject {
    double Q2, nu, Xb, y, W;
    ClassDef(DISVars, 1);
};
struct EPPI0Vars : public TObject {
    int    detPi0, numProtons, numElectrons, numPhotons;
    double pi0_p, pi0_theta, pi0_phi;
    double t, trentoPhi;
    double m_gg, m2_miss, m2_epX, m2_epi0X;
    double E_miss, px_miss, py_miss, pz_miss, pT_miss;
    ClassDef(EPPI0Vars, 1);
};
