#include "Kinematics.h"
#include <cmath>
#include <algorithm> // std::max

Kinematics::Kinematics(const TLorentzVector& e_f, double E_beam, double m_target) : 

    e_f_(e_f), E_beam_(E_beam), m_target_(m_target) {

        //std::cout << "Kinematics received e_f = " <<  e_f.Px() << ", " << e_f.Py() << ", " << e_f.Pz() << ", " << e_f.E() <<
        //                                            ", E_beam = " << E_beam << ", m_target = " << m_target << std::endl;
        beam_.SetPxPyPzE(0, 0, E_beam, E_beam);
        target_.SetPxPyPzE(0, 0, 0, m_target);
        computeDIS();
        //std::cout << "Kinematics computeDIS() obtained Q2 = " << Q2() << ", Xb = " << Xb() << std::endl;
    }

Kinematics::Kinematics(const TLorentzVector& e_f,
    const TLorentzVector& p_f,
    const TLorentzVector& pi0,
    double E_beam,
    double m_target) :

    e_f_(e_f), p_f_(p_f), pi0_(pi0), E_beam_(E_beam), m_target_(m_target) {

        beam_.SetPxPyPzE(0, 0, E_beam, E_beam);
        target_.SetPxPyPzE(0, 0, 0, m_target);
        computeDIS();
        computeExclusivePi0();
    }

void Kinematics::computeDIS() {
    TLorentzVector q = beam_ - e_f_;
    Q2_ = -q.Mag2();                    // Q² = -q²
    nu_ = q.E();                        // Energy transfer
    Xb_ = -q.Mag2() / (2 * m_target_ * q.E());   // Bjorken x
    y_  = q.E() / beam_.E(); // fractional energy loss via virtual photon
    W_  = sqrt(std::max(0.0, m_target_ * m_target_ + 2 * m_target_ * q.E() + q.Mag2()));   
}

void Kinematics::computeExclusivePi0() {
    // Missing momentum: beam + target - e' - p - pi0
    missing_       = beam_ + target_ - e_f_ - p_f_ - pi0_;
    missing_epX_   = beam_ + target_ - e_f_ - p_f_;
    missing_epi0X_ = beam_ + target_ - e_f_ - pi0_;

    // Mandelstam t = (p' - p)^2 = (final proton - target)^2. NOTE: t is redundant and already computed when proton is detected.
    TLorentzVector t = p_f_ - target_;
    t_ = t.Mag2();  

    // Trento Phi: angle from lepton plane to hadron plane, w.r.t. virtual photon momentum
    // 3 vectors
    TVector3 e_i = beam_.Vect();     // incoming electron
    TVector3 e_f = e_f_.Vect(); // scattered electron
    TVector3 q   = (e_i - e_f);     // virtual photon
    TVector3 h_f = p_f_.Vect();  // final hadron, e.g., proton

    // Plane normals (order matters!)
    TVector3 n_lepton = e_i.Cross(e_f).Unit();     // lepton plane
    TVector3 n_hadron = h_f.Cross(q).Unit();        // hadron plane

    double cosPhi = n_lepton.Dot(n_hadron);
    double sinPhi = q.Unit().Dot(n_lepton.Cross(n_hadron));
    phiT_ = std::atan2(sinPhi, cosPhi);

}
    