#include "BranchVars.h"

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

template <typename T>
inline double safeGetIfEnabled(const std::unordered_set<std::string>& enabledBranches, const std::string& branch, T expr) {
    if (!enabledBranches.count(branch)) return -999;
    double val = static_cast<double>(expr);
    return !std::isnan(val) ? val : NAN;
}

template <typename T>
inline double safeGetIfDetIfEnabled(int det, int expectedDet, const std::unordered_set<std::string>& enabledBranches, const std::string& branch, T expr) {
    if (!enabledBranches.count(branch)) return -999;
    double val = static_cast<double>(expr);
    return (det == expectedDet && !std::isnan(val)) ? val : NAN;
}

// Based on particle status, returns int representation of FT (0), FD (1), or CD (2)
int getDetector(int status) {
    const int absStatus = std::abs(status);

    if (absStatus >= 1000 && absStatus < 2000) return 0; // FT
    if (absStatus >= 2000 && absStatus < 4000) return 1; // FD
    if (absStatus >= 4000 && absStatus < 5000) return 2; // CD

    return -999; // Unknown detector
}

void EventVars::fill(const std::unordered_set<std::string>& enabledBranches, clas12::clas12reader& c12) {
    runNum       = c12.runconfig()->getRun();
    eventNum     = c12.runconfig()->getEvent();
    helicity     = c12.event()->getHelicity();
    charge       = c12.event()->getBeamCharge();

    numElectrons = safeGetIfEnabled(enabledBranches, "numElectrons", c12.getByID(11).size());
    numProtons   = safeGetIfEnabled(enabledBranches, "numProtons", c12.getByID(2212).size());
    numPhotons   = safeGetIfEnabled(enabledBranches, "numPhotons", c12.getByID(22).size());

}

void RecVars::fill(const std::unordered_set<std::string>& enabledBranches, clas12::region_particle* rec) {
    pid           = rec->getPid();
    charge        = rec->par()->getCharge();
    status        = rec->par()->getStatus();
    det           = getDetector(rec->par()->getStatus());
    sector        = rec->getSector();
    p             = rec->getP();
    beta          = rec->par()->getBeta();
    theta         = rec->getTheta();
    phi           = rec->getPhi();
    chi2pid       = rec->par()->getChi2Pid();

    px            = safeGetIfEnabled(enabledBranches, "px", rec->par()->getPx());
    py            = safeGetIfEnabled(enabledBranches, "py", rec->par()->getPy());
    pz            = safeGetIfEnabled(enabledBranches, "pz", rec->par()->getPz());

    vx            = safeGetIfEnabled(enabledBranches, "vx", rec->par()->getVx());
    vy            = safeGetIfEnabled(enabledBranches, "vy", rec->par()->getVy());
    vz            = safeGetIfEnabled(enabledBranches, "vz", rec->par()->getVz());
    time          = safeGetIfEnabled(enabledBranches, "time", rec->getTime());

    xFT           = safeGetIfDetIfEnabled(det, 0, enabledBranches, "xFT", rec->ft(FTCAL)->getX());
    yFT           = safeGetIfDetIfEnabled(det, 0, enabledBranches, "yFT", rec->ft(FTCAL)->getY());

    xDC1          = safeGetIfDetIfEnabled(det, 1, enabledBranches, "xDC1", rec->traj(DC, 6)->getX());
    yDC1          = safeGetIfDetIfEnabled(det, 1, enabledBranches, "yDC1", rec->traj(DC, 6)->getY());
    xDC2          = safeGetIfDetIfEnabled(det, 1, enabledBranches, "xDC2", rec->traj(DC, 18)->getX());
    yDC2          = safeGetIfDetIfEnabled(det, 1, enabledBranches, "yDC2", rec->traj(DC, 18)->getY());
    xDC3          = safeGetIfDetIfEnabled(det, 1, enabledBranches, "xDC3", rec->traj(DC, 36)->getX());
    yDC3          = safeGetIfDetIfEnabled(det, 1, enabledBranches, "yDC3", rec->traj(DC, 36)->getY());

    uPCAL         = safeGetIfDetIfEnabled(det, 1, enabledBranches, "uPCAL", rec->cal(1)->getLu());
    vPCAL         = safeGetIfDetIfEnabled(det, 1, enabledBranches, "vPCAL", rec->cal(1)->getLv());
    wPCAL         = safeGetIfDetIfEnabled(det, 1, enabledBranches, "wPCAL", rec->cal(1)->getLw());

    uECIN         = safeGetIfDetIfEnabled(det, 1, enabledBranches, "uECIN", rec->cal(4)->getLu());
    vECIN         = safeGetIfDetIfEnabled(det, 1, enabledBranches, "vECIN", rec->cal(4)->getLv());
    wECIN         = safeGetIfDetIfEnabled(det, 1, enabledBranches, "wECIN", rec->cal(4)->getLw());

    uECOUT        = safeGetIfDetIfEnabled(det, 1, enabledBranches, "uECOUT", rec->cal(7)->getLu());
    vECOUT        = safeGetIfDetIfEnabled(det, 1, enabledBranches, "vECOUT", rec->cal(7)->getLv());
    wECOUT        = safeGetIfDetIfEnabled(det, 1, enabledBranches, "wECOUT", rec->cal(7)->getLw());

    E_PCAL        = safeGetIfDetIfEnabled(det, 1, enabledBranches, "E_PCAL",  rec->cal(1)->getEnergy());
    E_ECIN        = safeGetIfDetIfEnabled(det, 1, enabledBranches, "E_ECIN",  rec->cal(4)->getEnergy());
    E_ECOUT       = safeGetIfDetIfEnabled(det, 1, enabledBranches, "E_ECOUT", rec->cal(7)->getEnergy());

    edge_cvt1     = safeGetIfDetIfEnabled(det, 2, enabledBranches, "edge_cvt1",  rec->traj(CVT, 1)->getEdge());
    edge_cvt3     = safeGetIfDetIfEnabled(det, 2, enabledBranches, "edge_cvt3",  rec->traj(CVT, 3)->getEdge());
    edge_cvt5     = safeGetIfDetIfEnabled(det, 2, enabledBranches, "edge_cvt5",  rec->traj(CVT, 5)->getEdge());
    edge_cvt7     = safeGetIfDetIfEnabled(det, 2, enabledBranches, "edge_cvt7",  rec->traj(CVT, 7)->getEdge());
    edge_cvt12    = safeGetIfDetIfEnabled(det, 2, enabledBranches, "edge_cvt12", rec->traj(CVT, 12)->getEdge());

    double x_cvt  = safeGetIfDetIfEnabled(det, 2, enabledBranches, "theta_cvt", rec->traj(CVT, 1)->getX());
    double y_cvt  = safeGetIfDetIfEnabled(det, 2, enabledBranches, "theta_cvt", rec->traj(CVT, 1)->getY());
    double z_cvt  = safeGetIfDetIfEnabled(det, 2, enabledBranches, "theta_cvt", rec->traj(CVT, 1)->getZ());

    double r_cvt  = (!std::isnan(x_cvt) && !std::isnan(y_cvt) && !std::isnan(z_cvt))
                        ? std::sqrt(x_cvt * x_cvt + y_cvt * y_cvt + z_cvt * z_cvt) : NAN;

    theta_cvt = (!std::isnan(r_cvt) && enabledBranches.count("theta_cvt")) 
                    ? std::acos(z_cvt / r_cvt) : -999;

    phi_cvt = (!std::isnan(x_cvt) && !std::isnan(y_cvt) && enabledBranches.count("phi_cvt")) 
                    ? std::atan2(y_cvt, x_cvt) : -999;
}

void RecVars::fill(const std::unordered_set<std::string>& enabledBranches, clas12::region_particle* rec,
                   double pIn, double thetaIn, double phiIn) { 
    fill(enabledBranches, rec);

    // Apply corrections
    p     = pIn;
    theta = thetaIn;
    phi   = phiIn;

    px = safeGetIfEnabled(enabledBranches, "px", p * std::sin(theta) * std::cos(phi));
    py = safeGetIfEnabled(enabledBranches, "py", p * std::sin(theta) * std::sin(phi));
    pz = safeGetIfEnabled(enabledBranches, "pz", p * std::cos(theta));
    
}

void GenVars::fill(clas12::mcparticle* mcParticles, int j) {
    pid  = mcParticles->getPid(j);
    TVector3 v_gen = TVector3(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j));
    p = v_gen.Mag();
    theta = v_gen.Theta();
    phi = v_gen.Phi();
}

void DISVars::fill(const TLorentzVector& lv_ePrime, double ebeam, double m_target) {
    // Beam and target 4-vectors:
    TLorentzVector lv_beam(0, 0, ebeam, ebeam);
    TLorentzVector lv_target(0, 0, 0, m_target);
    
    // Four-momentum transfer q
    TLorentzVector q = lv_beam - lv_ePrime;

    // DIS variables
    Q2  = -q.M2();                             // GeV^2
    nu  = ebeam - lv_ePrime.E();               // Energy transfer (GeV)
    Xb  = Q2 / (2.0 * m_target * nu);          // Bjorken x
    y   = nu / ebeam;                          // Inelasticity
    double W2  = (lv_target + q).M2();            // Invariant mass squared of hadronic system
    W   = std::sqrt(std::max(0.0, W2));

};

double EPPI0Vars::computeTrentoPhi(const TLorentzVector& lv_beam,
                                   const TLorentzVector& lv_target,
                                   const TLorentzVector& lv_ePrime,
                                   const TLorentzVector& lv_hadron) const {
    // Virtual photon
    TLorentzVector lv_q = lv_beam - lv_ePrime;

    // Total gamma*-p system
    TLorentzVector lv_tot = lv_q + lv_target;

    // Boost vector to gamma*-p CM
    TVector3 beta = -lv_tot.BoostVector();

    // Boost all relevant vectors
    TLorentzVector beam_b   = lv_beam;   beam_b.Boost(beta);
    TLorentzVector eout_b   = lv_ePrime; eout_b.Boost(beta);
    TLorentzVector q_b      = lv_q;      q_b.Boost(beta);
    TLorentzVector had_b    = lv_hadron; had_b.Boost(beta);

    // Compute plane normals
    TVector3 n_lepton = beam_b.Vect().Cross(eout_b.Vect());
    if(n_lepton.Mag2() > 0) n_lepton = n_lepton.Unit();

    TVector3 n_hadron = had_b.Vect().Cross(q_b.Vect());
    if(n_hadron.Mag2() > 0) n_hadron = n_hadron.Unit();

    // Signed Trento phi
    TVector3 q_hat = q_b.Vect().Unit();
    double cosPhi = n_lepton.Dot(n_hadron);
    double sinPhi = q_hat.Dot(n_lepton.Cross(n_hadron));

    double phi = std::atan2(sinPhi, cosPhi); // [-pi, pi]

    return phi;
}


void EPPI0Vars::fill(const TLorentzVector& lv_ePrime, const TLorentzVector& lv_pPrime, 
                     const TLorentzVector& lv_obj1, const TLorentzVector& lv_obj2, double ebeam, bool is_rad, double m_target) {
    // Beam and target 4-vectors:
    TLorentzVector lv_beam(0, 0, ebeam, ebeam);
    TLorentzVector lv_target(0, 0, 0, m_target);
    TLorentzVector lv_pi0;
    TLorentzVector lv_g1, lv_g2;        // always stored (rad → zero)
    TLorentzVector lv_g_rad;            // always stored (non-rad → zero)

    if (!is_rad) {
        // ---------------------------
        // NON-RADIATIVE: obj1 = g1, obj2 = g2
        // ---------------------------
        lv_g1  = lv_obj1;
        lv_g2  = lv_obj2;
        lv_pi0 = lv_obj1 + lv_obj2;

        // No radiative photon
        lv_g_rad.SetPxPyPzE(0,0,0,0);
    }
    else {
        // ---------------------------
        // RADIATIVE: obj1 = pi0, obj2 = gamma_rad
        // ---------------------------
        lv_pi0   = lv_obj1;   // generator π0
        lv_g_rad = lv_obj2;   // radiative photon

        // No decay photons stored
        lv_g1.SetPxPyPzE(0,0,0,0);
        lv_g2.SetPxPyPzE(0,0,0,0);
    }


    TLorentzVector lv_missing = lv_beam + lv_target - lv_ePrime - lv_pPrime - lv_pi0;
    TLorentzVector lv_epX     = lv_beam + lv_target - lv_ePrime - lv_pPrime;
    TLorentzVector lv_epi0X   = lv_beam + lv_target - lv_ePrime - lv_pi0;

    pi0_p         = lv_pi0.P();
    pi0_theta     = lv_pi0.Theta();
    pi0_phi       = lv_pi0.Phi();
    pi0_deltaPhi  = TVector2::Phi_mpi_pi(lv_pi0.Phi() - lv_epX.Phi());
    pi0_thetaX    = lv_pi0.Vect().Angle(lv_epX.Vect());
    m_gg          = lv_pi0.M();

    m2_miss       = lv_missing.M2();
    m2_epX        = lv_epX.M2();
    m2_epi0X      = lv_epi0X.M2();
    m_eggX        = (m2_epi0X >= 0.0) ? std::sqrt(m2_epi0X) : NAN;
    E_miss        = lv_missing.E();
    px_miss       = lv_missing.Px();
    py_miss       = lv_missing.Py();
    pz_miss       = lv_missing.Pz();
    pT_miss       = lv_missing.Pt();

    theta_e_g1  = is_rad ? NAN : lv_ePrime.Vect().Angle(lv_g1.Vect());
    theta_e_g2  = is_rad ? NAN : lv_ePrime.Vect().Angle(lv_g2.Vect());
    theta_g1_g2 = is_rad ? NAN : lv_g1.Vect().Angle(lv_g2.Vect());

    // Radiative photon quantities
    g_rad_E     = is_rad ? lv_g_rad.E() : NAN;
    g_rad_p     = is_rad ? lv_g_rad.P() : NAN;
    g_rad_theta = is_rad ? lv_g_rad.Theta() : NAN;
    g_rad_phi   = is_rad ? lv_g_rad.Phi() : NAN;

    // NOTE: "t" is implicitly -t!
    t = -1 * (lv_target - lv_pPrime).M2();

    // // Trento Phi: angle from lepton plane to hadron plane, w.r.t. virtual photon momentum
    // // 3 vectors
    // TVector3 v_beam      = lv_beam.Vect();            // incoming electron
    // TVector3 v_ePrime    = lv_ePrime.Vect();         // scattered electron
    // TVector3 v_q         = (v_beam - v_ePrime);     // virtual photon
    // TVector3 v_pPrime    = lv_pPrime.Vect();         // final hadron, e.g., proton

    // // Plane normals (order matters!)
    // TVector3 n_lepton = v_beam.Cross(v_ePrime).Unit();     // lepton plane
    // TVector3 n_hadron = v_pPrime.Cross(v_q).Unit();        // hadron plane

    // double cosPhi = n_lepton.Dot(n_hadron);
    // double sinPhi = v_q.Unit().Dot(n_lepton.Cross(n_hadron));

    // trentoPhi = std::atan2(sinPhi, cosPhi);

    // Using recoil proton
    trentoPhi = computeTrentoPhi(lv_beam, lv_target, lv_ePrime, lv_pPrime);

};

// Implement ClassImp macros for each class to generate ROOT dictionary info
ClassImp(EventVars);
ClassImp(RecVars);
ClassImp(GenVars);
ClassImp(DISVars);
ClassImp(EPPI0Vars);