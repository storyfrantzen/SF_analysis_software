#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "TLorentzVector.h"
#include "PhysicalConstants.h"

class Kinematics {
public:
    // DIS-only constructor which takes final electron `e_f`, default (6.535 GeV) beam, and default (proton) target
    Kinematics(const TLorentzVector& e_f, double ebeam = 6.535, double m_target = PROTON_MASS);

    // EPPI0-specific constructor for DIS + exclusivity analysis
    Kinematics(const TLorentzVector& e_f,
               const TLorentzVector& p_f,
               const TLorentzVector& pi0,
               double ebeam = 6.535,
               double m_target = PROTON_MASS);

    // Accessors for DIS kinematics:
    double Q2() const { return Q2_; }
    double nu() const { return nu_; }
    double Xb() const { return Xb_; }
    double y()  const { return y_;  }
    double W()  const { return W_;  }

    // Exclusivity accessors:

    bool channelCheck() const { return (Q2_ >= 1.0 && W_ >= 2.0 && y_ <= 0.80); }

    double t() const { return t_; }
    double missingM() const { return missing_.M(); }
    double phiT() const { return phiT_; }

    TLorentzVector missingP4() const { return missing_; }  // missing momentum expressed as 4-vector
    TLorentzVector lv_epX()    const { return missing_epX_; } // final state has e' and p' but misses (pi0)', i.e., ep -> e'p'X
    TLorentzVector lv_epi0X()  const { return missing_epi0X_; } // ep -> e'(pi0)'X
    
    // Access lv_beam and lv_target:
    TLorentzVector Beam()   const { return beam_; }
    TLorentzVector Target() const { return target_; }

private:

    double ebeam_, m_target_;

    // Possible 4-vectors which appear in an event: 
    TLorentzVector beam_, target_, e_f_, p_f_, pi0_;

    // DIS quantities:
    double Q2_, nu_, Xb_, y_, W_;
    void computeDIS();  // Internal method to update Q2, Nu, Xb, W

    // Exclusivity quantities:
    TLorentzVector missing_, missing_epX_, missing_epi0X_;
    double t_;
    double phiT_;
    void computeExclusivePi0();
    
};

#endif
