#include "Vars.h"

// template <typename T>
// void registerBranchesWithPrefix(std::vector<BranchInfo>& list, const T& obj, const std::string& prefix) {
//     // [Not in use] helper function to streamline adding arbitrarily many particles to the supported Master Branches List
//     registerBranchInfo(list, prefix + "pid",    &const_cast<T&>(obj).pid,    "I");
//     registerBranchInfo(list, prefix + "charge", &const_cast<T&>(obj).charge, "I");
//     // repeat for all members...
// }

template <typename T>
inline double safeGet(T expr) {
    double val = static_cast<double>(expr);
    return !std::isnan(val) ? val : NAN;
}

template <typename T>
inline double safeGetIfDet(int det, int expectedDet, T expr) {
    double val = static_cast<double>(expr);
    return (det == expectedDet && !std::isnan(val)) ? val : NAN;
}

// Implement ClassImp macros for each class to generate ROOT dictionary info
ClassImp(EventVars);
ClassImp(RecVars);
ClassImp(GenVars);
ClassImp(DISVars);
ClassImp(EPPI0Vars);

// Based on particle status, returns int representation of FT (0), FD (1), or CD (2)
int getDetector(int status) {
    const int absStatus = std::abs(status);

    if (absStatus >= 1000 && absStatus < 2000) return 0; // FT
    if (absStatus >= 2000 && absStatus < 4000) return 1; // FD
    if (absStatus >= 4000 && absStatus < 5000) return 2; // CD

    return -999; // Unknown detector
}

void EventVars::fill(clas12::clas12reader& c12) {
        runNum = c12.runconfig()->getRun();
        eventNum = c12.runconfig()->getEvent();
        helicity = c12.event()->getHelicity();
    }

void RecVars::fill(clas12::region_particle* part) {
    pid             = part->getPid();
    charge          = part->par()->getCharge();
    status          = part->par()->getStatus();
    det             = getDetector(part->par()->getStatus());
    sector          = part->getSector();

    p               = part->getP();
    px              = part->par()->getPx();
    py              = part->par()->getPy();
    pz              = part->par()->getPz();
    beta            = part->par()->getBeta();
    theta           = part->getTheta();
    phi             = part->getPhi();
    vx              = part->par()->getVx();
    vy              = part->par()->getVy();
    vz              = part->par()->getVz();
    chi2pid         = part->par()->getChi2Pid();
    time            = part->getTime();

    xFT             = safeGetIfDet(det, 0, part->ft(FTCAL)->getX());
    yFT             = safeGetIfDet(det, 0, part->ft(FTCAL)->getY());
    xDC1            = safeGetIfDet(det, 1, part->traj(DC, 6)->getX());
    yDC1            = safeGetIfDet(det, 1, part->traj(DC, 6)->getY());
    xDC2            = safeGetIfDet(det, 1, part->traj(DC, 18)->getX());
    yDC2            = safeGetIfDet(det, 1, part->traj(DC, 18)->getY());
    xDC3            = safeGetIfDet(det, 1, part->traj(DC, 36)->getX());
    yDC3            = safeGetIfDet(det, 1, part->traj(DC, 36)->getY());

    uPCAL           = safeGetIfDet(det, 1, part->cal(1)->getLu());
    vPCAL           = safeGetIfDet(det, 1, part->cal(1)->getLv());
    wPCAL           = safeGetIfDet(det, 1, part->cal(1)->getLw());

    uECIN           = safeGetIfDet(det, 1, part->cal(4)->getLu());
    vECIN           = safeGetIfDet(det, 1, part->cal(4)->getLv());
    wECIN           = safeGetIfDet(det, 1, part->cal(4)->getLw());

    uECOUT          = safeGetIfDet(det, 1, part->cal(7)->getLu());
    vECOUT          = safeGetIfDet(det, 1, part->cal(7)->getLv());
    wECOUT          = safeGetIfDet(det, 1, part->cal(7)->getLw());

    E_PCAL          = safeGetIfDet(det, 1, part->cal(1)->getEnergy());
    E_ECIN          = safeGetIfDet(det, 1, part->cal(4)->getEnergy());
    E_ECOUT         = safeGetIfDet(det, 1, part->cal(7)->getEnergy());

    edge_cvt1       = safeGetIfDet(det, 2, part->traj(CVT, 1)->getEdge());
    edge_cvt3       = safeGetIfDet(det, 2, part->traj(CVT, 3)->getEdge());
    edge_cvt5       = safeGetIfDet(det, 2, part->traj(CVT, 5)->getEdge());
    edge_cvt7       = safeGetIfDet(det, 2, part->traj(CVT, 7)->getEdge());
    edge_cvt12      = safeGetIfDet(det, 2, part->traj(CVT, 12)->getEdge());

    double x_cvt    = safeGetIfDet(det, 2, part->traj(CVT, 1)->getX());
    double y_cvt    = safeGetIfDet(det, 2, part->traj(CVT, 1)->getY());
    double z_cvt    = safeGetIfDet(det, 2, part->traj(CVT, 1)->getZ());

    double r_cvt = (!std::isnan(x_cvt) && !std::isnan(y_cvt) && !std::isnan(z_cvt))
                    ? std::sqrt(x_cvt * x_cvt + y_cvt * y_cvt + z_cvt * z_cvt) : NAN;

    theta_cvt = (!std::isnan(r_cvt)) 
                    ? safeGetIfDet(det, 2, std::acos(z_cvt / r_cvt)) : NAN;

    phi_cvt = (!std::isnan(x_cvt) && !std::isnan(y_cvt))
                    ? safeGetIfDet(det, 2, std::atan2(y_cvt, x_cvt)) : NAN;
}
