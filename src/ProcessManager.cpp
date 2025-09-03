#include "ProcessManager.h"

#include "nlohmann/json.hpp"
#include "FiducialCuts.h"       // ProcessManager has a member instance of FC during filtering
#include "PhysicalConstants.h"  // contains useful constants
#include "clas12reader.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace clas12;

ProcessManager::ProcessManager(const nlohmann::json& config) { 

    // --- Basic Config -- 
    ebeam_       = config["ebeam"];
    torus_       = config["torus"];
    channel_     = config["channel"];
    outPrefix_   = config.value("outPrefix", "");

    // --- Fiducial Cuts ---
    FC_ = std::make_unique<FiducialCuts>();
    for (const auto& cut : config.value("fiducialCuts", std::vector<std::string>{})) { FC_->addCut(cut); }

    // --- Kinematic Corrections ---    
    if (config.contains("kinCorrections")) {
        if (config["kinCorrections"].is_string()) {
            KC_ = std::make_unique<KinematicCorrections>(config["kinCorrections"].get<std::string>());
        } else if (config["kinCorrections"].is_object()) {
            KC_ = std::make_unique<KinematicCorrections>(config["kinCorrections"]);
        } else {
            KC_ = std::make_unique<KinematicCorrections>(nlohmann::json{}); 
        }
    } else KC_ = std::make_unique<KinematicCorrections>(nlohmann::json{});

    // --- SF Cuts ---
    if (config.contains("SFCuts")) {
        if (config["SFCuts"].is_string()) {
            SF_ = std::make_unique<SFCuts>(config["SFCuts"].get<std::string>());
        } else if (config["SFCuts"].is_object()) {
            SF_ = std::make_unique<SFCuts>(config["SFCuts"]);
        } else {
            SF_ = std::make_unique<SFCuts>(nlohmann::json{}); // empty → no cuts
        }
    } else SF_ = std::make_unique<SFCuts>(nlohmann::json{}); // default pass-through
    

    std::stringstream ss;
    ss << "output/";
    if (!outPrefix_.empty()) ss << outPrefix_ << "_";
    ss << currentTimestamp() << "_" << ebeam_ << "_tor" << torus_ << "_" 
       << channel_ << ".root";
    std::string outFileName = ss.str();
    
    outFile_ = new TFile(outFileName.c_str(), "RECREATE");
    if (!outFile_ || outFile_->IsZombie()) {
        std::cerr << "Error opening output file: " << outFileName << std::endl;
        outFile_ = nullptr;
        return;
    }
    outFile_->cd();

    tree_ = new TTree("Events", "CLAS12 event data");
    tree_->SetAutoSave(500000);  // Autosave every 500k bytes or choose a suitable value
    tree_->SetAutoFlush(5000);  

    tree_->Branch("event", "EventVars", &ev_);
    tree_->Branch("dis", "DISVars", &dis_);

    if (channel_ == "genOnly") {
        tree_->Branch("gen_dis", "GenVars", &gen_dis_);
        tree_->Branch("gen", "GenVars", &gen_);
    } else if (channel_ == "inclusiveRec") {
        tree_->Branch("rec", "RecVars", &rec_);
    } else if (channel_ == "genMatch") {
        tree_->Branch("gen_dis",  "DISVars",     &gen_dis_);
        tree_->Branch("gen",      "GenVars",     &gen_);
        tree_->Branch("rec",      "RecVars",     &rec_);
    } else if (channel_ == "eppi0") {
        tree_->Branch("e",     "RecVars",     &e_);
        tree_->Branch("p",     "RecVars",     &p_);
        tree_->Branch("g",     "RecVars",     &g_);
        tree_->Branch("eppi0", "EPPI0Vars",   &eppi0_);
    }

    if (config.contains("branches") && !config["branches"].is_null()) {
        auto& b = config["branches"];
        for (const auto& var : b.value("event", std::vector<std::string>{})) enabledEvBranches_.insert(var);
        for (const auto& var : b.value("rec", std::vector<std::string>{})) enabledRecBranches_.insert(var);
        for (const auto& var : b.value("electron", std::vector<std::string>{})) enabledEleBranches_.insert(var);
        for (const auto& var : b.value("proton", std::vector<std::string>{})) enabledProBranches_.insert(var);
        for (const auto& var : b.value("photon", std::vector<std::string>{})) enabledPhoBranches_.insert(var);
    } 
}

//////////////////////////////////////////////////////////
// Helper functions called by bigger functions: ///
//////////////////////////////////////////////////////////

std::string ProcessManager::currentTimestamp() const {
    auto now = std::chrono::system_clock::now();
    std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&tt), "%m%d_%H%M");
    return ss.str();
}

// Based on particle status, returns int representation of FT (0), FD (1), or CD (2)
int ProcessManager::getDetector(int status) {
        const int absStatus = std::abs(status);

        if (absStatus >= 1000 && absStatus < 2000) return 0; // FT
        if (absStatus >= 2000 && absStatus < 4000) return 1; // FD
        if (absStatus >= 4000 && absStatus < 5000) return 2; // CD

        return -999; // Unknown detector
}

bool ProcessManager::channelCheck(float Q2, float W, float y) {
    return Q2 >= 1 && W >= 2 && y <= 0.8;
}

// Currently in use. Default is -8 <= z <= 2.
bool ProcessManager::passesVertexCut(const float vz, const int zmin, const int zmax) { 
    return vz >= zmin && vz <= zmax; 
}

void ProcessManager::finalize() {
    if (outFile_) {
        outFile_->cd();
        if (tree_) {
            auto bytesWritten = tree_->Write();
            if (bytesWritten > 0) {
                std::cout << "[ProcessManager] Successfully saved tree to file: "
                          << outFile_->GetName()
                          << " (" << bytesWritten << " bytes written)" << std::endl;
            } else {
                std::cerr << "[ProcessManager] Error: Failed to write tree to file: "
                          << outFile_->GetName() << std::endl;
            }
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

auto getPhotonCalE = [](clas12::region_particle* g) {
    int det = ProcessManager::getDetector(g->par()->getStatus());
    if (det == 0) {
        return g->ft(FTCAL) ? g->ft(FTCAL)->getEnergy() : 0.0;
    } else if (det == 2) {
        return 0.0;
    }
    return (g->cal(1) ? g->cal(1)->getEnergy() : 0.0) + 
           (g->cal(4) ? g->cal(4)->getEnergy() : 0.0) + 
           (g->cal(7) ? g->cal(7)->getEnergy() : 0.0);
};

//////////////////////////////////////////////////////////
// MAIN WORKFLOW PERFORMED HERE: //
//////////////////////////////////////////////////////////

void ProcessManager::processEvent(clas12::clas12reader& c12) {
    eventsProcessed_++;
    if (channel_ == "genOnly") {

        // Generated particles
        auto mcParticles = c12.mcparts();
        int numGen = mcParticles->getRows();

        bool electronFound = false;  
        TLorentzVector lv_ePrime; // will store the first electron’s 4-vector

        for (int j = 0; j < numGen; j++) {
            // initialize ALL vars:
            ev_.flush();
            dis_.flush();
            gen_.flush();

            ev_.fill(enabledEvBranches_, c12);
            gen_.fill(mcParticles, j);

            int pid = mcParticles->getPid(j);
            TVector3 v_p(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j));

            if (pid == 11 && !electronFound) {
                // Calculate and store only for the first electron
                lv_ePrime.SetVectM(v_p, ELECTRON_MASS);
                electronFound = true;
            } else if (!electronFound) {
                dis_.flush();
            } else dis_.fill(lv_ePrime, ebeam_);

            tree_->Fill();
            numFills_++;
        }

    } else if (channel_ == "inclusiveRec") {
        
        auto recParticles = c12.getDetParticles();   
        int numRec = recParticles.size();     

        bool electronFound = false;  
        TLorentzVector lv_ePrime; // will store the first electron’s 4-vector       

        for (int i = 0; i < numRec; i++) {
            // Initialize ALL vars:
            ev_.flush();
            dis_.flush();
            rec_.flush();

            // Begin:
            auto& p = recParticles[i];

            int pid  = p->getPid();
            int det  = getDetector(p->par()->getStatus());
            int sector = p->getSector();
            float vz = p->par()->getVz();

            if (pid == 11) {
                if (p->getP() < 1) continue;

                bool passesFC = ((det == 0 && FC_->passesFT(p)) || 
                                 (det == 1 && FC_->passesDC(p, torus_) && FC_->passesECAL(p)));
                if (!(passesFC && passesVertexCut(vz))) {
                    if (i == 0) return;
                    else continue;
                }
                
                if (det == 1) {
                    float ePCAL  = p->cal(1) ? p->cal(1)->getEnergy() : 0;
                    float eECIN  = p->cal(4) ? p->cal(4)->getEnergy() : 0;
                    float eECOUT = p->cal(7) ? p->cal(7)->getEnergy() : 0;
                    float sf = (ePCAL + eECIN + eECOUT) / p->getP();
                    if (!SF_->passTriangleCut(ePCAL, eECIN, p->getP())) continue;
                    if (!SF_->passSigmaCut(sector, sf, p->getP())) continue;
                }

                if (!electronFound) {
                    lv_ePrime.SetPxPyPzE(p->par()->getPx(), 
                                         p->par()->getPy(), 
                                         p->par()->getPz(), 
                                         std::sqrt(ELECTRON_MASS * ELECTRON_MASS + p->getP() * p->getP()));
                    electronFound = true;
                } 
            }
            else if (pid == 2212) {
                bool passesFC = ((det == 1 && FC_->passesDC(p, torus_)) || 
                                 (det == 2 && FC_->passesCVT(p)));
                if (!(passesFC && passesVertexCut(vz))) continue;
            }
            else if (pid == 22) {
                if (p->par()->getBeta() < 0.9 || p->par()->getBeta() > 1.1) continue;
                bool passesFC = (det == 0 && FC_->passesFT(p)) || (det == 1 && FC_->passesECAL(p));
                if (!passesFC) continue;
            }

            // If electron already found, repeat DIS fill for all subsequent particles
            if (electronFound) dis_.fill(lv_ePrime, ebeam_);
            ev_.fill(enabledEvBranches_, c12);
            if (pid == 2212) {
                double p_corr = p->getP() + KC_->deltaP(p->getP(), p->getTheta(), det == 1);
                double theta_corr = p->getTheta() + KC_->deltaTheta(p->getP(), p->getTheta(), det == 1);
                double phi_wrap = p->getPhi() < 0 ? p->getPhi() + 2*M_PI : p->getPhi();
                double phi_corr = phi_wrap + KC_->deltaPhi(p->getP(), p->getTheta(), det == 1);

                if (phi_corr > M_PI) phi_corr -= 2*M_PI;
                if (phi_corr < -M_PI) phi_corr += 2*M_PI;
                rec_.fill(enabledRecBranches_, p, p_corr, theta_corr, phi_corr);
            } else rec_.fill(enabledRecBranches_, p);
            tree_->Fill();
            numFills_++;
        }
    }
    else if (channel_ == "genMatch") {

        // Initialize ALL vars:
            ev_.flush();
            dis_.flush();
            gen_dis_.flush();
            rec_.flush();
            gen_.flush();
        
        auto mcParticles = c12.mcparts();
        auto recParticles = c12.getDetParticles();

        int numGen = mcParticles->getRows(); 
        int numRec = recParticles.size();      

        for (int j = 0; j < numGen; j++) {
            int gen_pid = mcParticles->getPid(j);
            TVector3 v_gen_p(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j));

            double bestDeltaTheta = 999.0;
            int bestMatchIndex = -1;

            for (int i = 0; i < numRec; i++) {
                auto& rec = recParticles[i];
                if (rec->par()->getPid() != gen_pid) continue;

                TVector3 v_rec_p(rec->par()->getPx(), rec->par()->getPy(), rec->par()->getPz());
                // double deltaPhi = fabs(v_gen_p.Phi() - v_rec_p.Phi()) * 180/M_PI;
                // double deltaTheta = fabs(v_gen_p.Theta() - v_rec_p.Theta()) * 180/M_PI;
                double deltaTheta = v_gen_p.Angle(v_rec_p) * 180/M_PI; // in degrees

                if (deltaTheta < bestDeltaTheta) {
                    bestDeltaTheta = deltaTheta;
                    bestMatchIndex = i;
                }
            }

            // Threshold for match quality:
            if (bestMatchIndex < 0 || bestDeltaTheta > 3) continue;

            auto& rec  = recParticles[bestMatchIndex];
            int det    = getDetector(rec->par()->getStatus());
            int sec    = rec->getSector();
            float vz   = rec->par()->getVz();

            bool electronFound = false;  
            TLorentzVector lv_gen_ePrime;  
            TLorentzVector lv_ePrime; // will store the first electron’s 4-vector       

            if (gen_pid == 11) {
                // Reject if electron has momentum less than threshold
                if (rec->getP() < 1) continue;

                bool passesFC = ((det == 0 && FC_->passesFT(rec)) || 
                                 (det == 1 && FC_->passesDC(rec, torus_) && FC_->passesECAL(rec)));
                if (!(passesFC && passesVertexCut(vz))) continue;

                if (det == 1) {
                    float ePCAL  = rec->cal(1) ? rec->cal(1)->getEnergy() : 0;
                    float eECIN  = rec->cal(4) ? rec->cal(4)->getEnergy() : 0;
                    float eECOUT = rec->cal(7) ? rec->cal(7)->getEnergy() : 0;
                    float sf = (ePCAL + eECIN + eECOUT) / rec->getP();
                    if (!SF_->passTriangleCut(ePCAL, eECIN, rec->getP())) continue;
                    if (!SF_->passSigmaCut(sec, sf, rec->getP())) continue;
                }

                if (!electronFound) {
                    lv_gen_ePrime.SetPxPyPzE(v_gen_p.Px(), 
                                             v_gen_p.Py(), 
                                             v_gen_p.Pz(), 
                                             std::sqrt(ELECTRON_MASS * ELECTRON_MASS + v_gen_p.Mag() * v_gen_p.Mag()));

                    lv_ePrime.SetPxPyPzE(rec->par()->getPx(), 
                                         rec->par()->getPy(), 
                                         rec->par()->getPz(), 
                                         std::sqrt(ELECTRON_MASS * ELECTRON_MASS + rec->getP() * rec->getP()));
                    electronFound = true;
                } 
            } 

            else if (gen_pid == 2212) {
                // Reject if proton has momentum less than threshold (Yijie + Bobby)
                if (rec->getP() < 0.3) continue;

                bool passesFC = ((det == 1 && FC_->passesDC(rec, torus_)) || 
                                 (det == 2 && FC_->passesCVT(rec)));
                if (!(passesFC && passesVertexCut(vz))) return;
            }

            else if (gen_pid == 22) {
                if (rec->par()->getBeta() < 0.9 || rec->par()->getBeta() > 1.1) continue;
                bool passesFC = (det == 0 && FC_->passesFT(rec)) || (det == 1 && FC_->passesECAL(rec));
                if (!passesFC) return;
            }

            // If electron already found, repeat DIS fill for all subsequent particles
            if (electronFound) {
                gen_dis_.fill(lv_gen_ePrime, ebeam_);
                dis_.fill(lv_ePrime, ebeam_);
            }

            ev_.fill(enabledEvBranches_, c12);
            gen_.fill(mcParticles, j);
            
            if (gen_pid == 2212) {
                double p_corr = rec->getP() + KC_->deltaP(rec->getP(), rec->getTheta(), det == 1);
                double theta_corr = rec->getTheta() + KC_->deltaTheta(rec->getP(), rec->getTheta(), det == 1);
                double phi_wrap = rec->getPhi() < 0 ? rec->getPhi() + 2*M_PI : rec->getPhi();
                double phi_corr = phi_wrap + KC_->deltaPhi(rec->getP(), rec->getTheta(), det == 1);

                if (phi_corr > M_PI) phi_corr -= 2*M_PI;
                if (phi_corr < -M_PI) phi_corr += 2*M_PI;
                rec_.fill(enabledRecBranches_, rec, p_corr, theta_corr, phi_corr);
            }
            else rec_.fill(enabledRecBranches_, rec);
            tree_->Fill();
            numFills_++;
        }
    }
    else if (channel_ == "eppi0") processEPPI0(c12);
}

void ProcessManager::processEPPI0(clas12::clas12reader& c12) {
    auto electrons = c12.getByID(11);   // PID 11 = electron
    auto protons   = c12.getByID(2212); // PID 2212 = proton
    auto photons   = c12.getByID(22);   // PID 22 = photon

    // Reject if not EPPI0 final state, according to eventbuilder:
    if (photons.size() < 2 || electrons.size() != 1 || protons.size() != 1) return;

    ev_.flush();
    e_.flush();
    p_.flush();
    g_.flush();
    dis_.flush();
    eppi0_.flush();

    clas12::region_particle* best_electron = nullptr;
    clas12::region_particle* best_proton   = nullptr;
    clas12::region_particle* best_g1       = nullptr;
    clas12::region_particle* best_g2       = nullptr;

    TLorentzVector lv_g1, lv_g2;
    int detPi0;

    double leastDeltaM = 999;

    for (const auto& ele : electrons) {

        // Reject if electron has momentum less than threshold
        if (ele->getP() < 1) continue;

        // Reject if electron vertex isn't in target fiducial volume
        if (!passesVertexCut(ele->par()->getVz())) continue;

        int detEle = getDetector(ele->par()->getStatus());
        int secEle = ele->getSector();
        // Reject if electron fails fiducial cuts (ONLY if cuts are listed!)
        if (detEle == 0 && !FC_->passesFT(ele)) continue;
        if (detEle == 1 && (!FC_->passesDC(ele, torus_) || !FC_->passesECAL(ele))) continue;

        if (detEle == 1) {
            float ePCAL  = ele->cal(1) ? ele->cal(1)->getEnergy() : 0;
            float eECIN  = ele->cal(4) ? ele->cal(4)->getEnergy() : 0;
            float eECOUT = ele->cal(7) ? ele->cal(7)->getEnergy() : 0;
            float sf = (ePCAL + eECIN + eECOUT) / ele->getP();
            if (!SF_->passTriangleCut(ePCAL, eECIN, ele->getP())) continue;
            if (!SF_->passSigmaCut(secEle, sf, ele->getP())) continue;
        }

        for (const auto& pro : protons) {

            // Reject if proton has momentum less than threshold (Yijie + Bobby)
            if (pro->getP() < 0.3) continue;

            // Reject if proton vertex isn't in target fiducial volume
            if (!passesVertexCut(pro->par()->getVz())) continue;

            int detPro = getDetector(pro->par()->getStatus());
            // Reject if proton fails fiducial cuts (ONLY if cuts are listed!)
            if (detPro == 1 && !FC_->passesDC(pro, torus_)) continue;
            if (detPro == 2 && !FC_->passesCVT(pro)) continue;

            for (size_t i = 0; i < photons.size(); ++i) {
                for (size_t j = i + 1; j < photons.size(); ++j) {

                    auto& g1 = photons[i];
                    auto& g2 = photons[j];

                    int det1 = getDetector(g1->par()->getStatus());
                    int det2 = getDetector(g2->par()->getStatus());

                    int sec1 = g1->getSector();
                    int sec2 = g2->getSector();

                    // Reject if either photon is detected in CD
                    if (det1 == 2 || det2 == 2) continue;

                    // Reject if photons are not detected in same sector
                    if (sec1 != sec2) continue;

                    // Reject if either photon has momentum less than threshold (Bobby's thesis):
                    if (g1->getP() < 0.4  || g2->getP() < 0.4) continue;

                    // Reject if either photon has unphysical beta
                    if (g1->par()->getBeta() < 0.9 || g1->par()->getBeta() > 1.1) continue;
                    if (g2->par()->getBeta() < 0.9 || g2->par()->getBeta() > 1.1) continue;

                    // Reject if either photon fails calorimeter energy threshold
                    if (getPhotonCalE(g1) < 0.15 || getPhotonCalE(g2) < 0.15) continue;

                    // Reject if either photon is detected in FT and fails FT fiducial cuts (fails ONLY if FTstandardCut is listed!)
                    if ((det1 == 0 && !FC_->passesFT(g1)) || (det2 == 0 && !FC_->passesFT(g2))) continue;

                    //Reject if either photon is detected in FD and fails ECAL fiducial cuts (fails ONLY if ECAL cuts are listed!)
                    if ((det1 == 1 && !FC_->passesECAL(g1)) || (det2 == 1 && !FC_->passesECAL(g2))) continue;

                    lv_g1.SetXYZM(g1->par()->getPx(), g1->par()->getPy(), g1->par()->getPz(), 0.0);
                    lv_g2.SetXYZM(g2->par()->getPx(), g2->par()->getPy(), g2->par()->getPz(), 0.0);
                    TLorentzVector lv_candidatePi0 = lv_g1 + lv_g2;

                    double deltaM = std::abs(lv_candidatePi0.M() - PI0_MASS);
                    if (deltaM > 0.05) continue;

                    if (deltaM < leastDeltaM) {
                        detPi0         = det1; // Assign to pi0 the detector of the photon
                        leastDeltaM    = deltaM;
                        best_electron  = ele;
                        best_proton    = pro;
                        best_g1        = g1;
                        best_g2        = g2;
                    }
                }
            }
        }
    }

    // Reject if no suitable EPPI0 candidate has been found, i.e., best_electron == nullptr:
    if (!best_electron) return; 

    ev_.fill(enabledEvBranches_, c12);

    // ELECTRON INFO:
    e_.fill(enabledEleBranches_, best_electron);

    // PROTON INFO:
    int detBestPro = getDetector(best_proton->par()->getStatus());
    double p_corr = best_proton->getP() + KC_->deltaP(best_proton->getP(), best_proton->getTheta(), detBestPro == 1);
    double theta_corr = best_proton->getTheta() + KC_->deltaTheta(best_proton->getP(), best_proton->getTheta(), detBestPro == 1);
    double phi_wrap = best_proton->getPhi() < 0 ? best_proton->getPhi() + 2*M_PI : best_proton->getPhi();
    double phi_corr = phi_wrap + KC_->deltaPhi(best_proton->getP(), best_proton->getTheta(), detBestPro == 1);
    if (phi_corr > M_PI) phi_corr -= 2*M_PI;
    if (phi_corr < -M_PI) phi_corr += 2*M_PI;
    p_.fill(enabledProBranches_, best_proton, p_corr, theta_corr, phi_corr);

    // PHOTON INFO:
    g_.fill(enabledPhoBranches_, best_g1);
    
    // DIS + EPPI0 INFO:
    
    // Electron 4-vector
    TLorentzVector lv_ePrime;
    lv_ePrime.SetPxPyPzE(best_electron->par()->getPx(),
                         best_electron->par()->getPy(),
                         best_electron->par()->getPz(),
                         std::sqrt(ELECTRON_MASS * ELECTRON_MASS + best_electron->getP() * best_electron->getP()));

    TLorentzVector lv_pPrime;
    lv_pPrime.SetPxPyPzE(best_proton->par()->getPx(),
                         best_proton->par()->getPy(),
                         best_proton->par()->getPz(),
                         std::sqrt(PROTON_MASS * PROTON_MASS + best_proton->getP() * best_proton->getP()));
    
    dis_.fill(lv_ePrime, ebeam_);

    lv_g1.SetXYZM(best_g1->par()->getPx(), best_g1->par()->getPy(), best_g1->par()->getPz(), 0.0);
    lv_g2.SetXYZM(best_g2->par()->getPx(), best_g2->par()->getPy(), best_g2->par()->getPz(), 0.0);

    eppi0_.fill(lv_ePrime, lv_pPrime, lv_g1, lv_g2, ebeam_);

    // EXCLUSIVITY CUTS:
    if (eppi0_.E_miss > 1) return;
    if (eppi0_.pT_miss > 0.2) return;
    if (eppi0_.theta_e_g1 * 180.0/M_PI < 4 || eppi0_.theta_e_g2 * 180.0/M_PI < 4) return;
    if (eppi0_.theta_g1_g2 * 180.0/M_PI < 1) return;

    // Reject if Q2 < 1, W < 2, or y > 0.8, i.e., not in standard DIS region
    if (!channelCheck(dis_.Q2, dis_.W, dis_.y)) return;

    numFills_++;

    // FILL TREE, ONCE PER EVENT:
    tree_->Fill();
}