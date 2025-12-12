#include "ProcessManager.h"
#include "nlohmann/json.hpp"
#include "FiducialCuts.h"       
#include "PhysicalConstants.h"  // useful constants
#include "clas12reader.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace clas12;

ProcessManager::ProcessManager(const nlohmann::json& config) { 

    // --- Basic Config -- //
    ebeam_       = config["ebeam"];
    torus_       = config["torus"];
    channel_     = config["channel"];
    tag_         = config.value("tag", "");

    // --- Fiducial Cuts --- //
    FC_ = std::make_unique<FiducialCuts>();
    for (const auto& cut : config.value("fiducialCuts", std::vector<std::string>{})) { FC_->addCut(cut); }

    // --- Kinematic Corrections --- //
    if (config.contains("kinCorrections")) {
        if (config["kinCorrections"].is_string()) {
            KC_ = std::make_unique<KinematicCorrections>(config["kinCorrections"].get<std::string>());
        } else if (config["kinCorrections"].is_object()) {
            KC_ = std::make_unique<KinematicCorrections>(config["kinCorrections"]);
        } else {
            KC_ = std::make_unique<KinematicCorrections>(nlohmann::json{}); 
        }
    } else KC_ = std::make_unique<KinematicCorrections>(nlohmann::json{});

    // --- SF Cuts --- //
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
    if (!tag_.empty()) ss << tag_ << "_";
    ss << channel_ << "_" << ebeam_ << "_tor" << torus_ << "_" << currentTimestamp() << ".root"; 
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

    if (channel_ == "inclusiveREC") {
        tree_->Branch("rec", "RecVars", &rec_);
    } else if (channel_ == "inclusiveMATCH") {
        tree_->Branch("gen_dis",  "DISVars",     &gen_dis_);
        tree_->Branch("rec",      "RecVars",     &rec_);
        tree_->Branch("gen",      "GenVars",     &gen_);
    } else if (channel_ == "eppi0REC") {
        tree_->Branch("e",     "RecVars",     &e_);
        tree_->Branch("p",     "RecVars",     &p_);
        tree_->Branch("g",     "RecVars",     &g_);
        tree_->Branch("eppi0", "EPPI0Vars",   &eppi0_);
    } else if (channel_ == "eppi0GEMC") {
        tree_->Branch("gen_dis", "DISVars", &gen_dis_);
        tree_->Branch("e",     "RecVars",     &e_);
        tree_->Branch("p",     "RecVars",     &p_);
        tree_->Branch("g",     "RecVars",     &g_);
        tree_->Branch("eppi0", "EPPI0Vars", &eppi0_);
        tree_->Branch("gen_eppi0", "EPPI0Vars", &gen_eppi0_);
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
// --- Helper functions called by bigger functions: --- //
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

bool ProcessManager::passesVertexCut(const int pid, const float vz, const float e_vz) { 
    if (ebeam_ > 10) {
        if (std::abs(pid) == 2212) {
            float diff = e_vz - vz;
            return std::abs(diff) < 20;
        }
        else if (pid == 11) {
            if (torus_ > 0) return vz >= -18 && vz <= 10;
            else return vz >= -13 && vz <= 12;
        }
        else return true;
    }
    else return vz >= -8 && vz <= 2; 
}

void ProcessManager::finalize(double totalCharge) {
    if (!outFile_) {
        std::cerr << "[ProcessManager] finalize(): No output file open!" << std::endl;
        return;
    }

    outFile_->cd();

    // 1. Write the main event tree if it exists
    if (tree_) {
        auto bytesWritten = tree_->Write();
        if (bytesWritten > 0) {
            std::cout << "[ProcessManager] Successfully saved main tree: "
                      << outFile_->GetName() << " (" << bytesWritten << " bytes written)" 
                      << std::endl;
        } else {
            std::cerr << "[ProcessManager] Error: Failed to write main tree!" << std::endl;
        }
        delete tree_;
        tree_ = nullptr;
    }

    // 2. Create a small summary tree with total charge and events processed
    TTree summary("Summary", "Summary information");
    summary.Branch("TotalCharge", &totalCharge, "TotalCharge/D");
    summary.Branch("EventsProcessed", &eventsProcessed_, "EventsProcessed/I");
    summary.Branch("Fills", &numFills_, "Fills/I");
    summary.Fill();
    summary.Write("", TObject::kOverwrite);

    std::cout << "[ProcessManager] Summary tree written: "
              << "TotalCharge = " << totalCharge << " nC, "
              << "EventsProcessed = " << eventsProcessed_ << ", "
              << "Fills = " << numFills_ << std::endl;

    // 3. Close the output file
    outFile_->Close();
    delete outFile_;
    outFile_ = nullptr;
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
//------------------- MAIN WORKFLOW: -------------------//
//////////////////////////////////////////////////////////

void ProcessManager::processEvent(clas12::clas12reader& c12) {
    eventsProcessed_++;
    bool fillTree = false;
    if (channel_ == "inclusiveREC") {
        
        auto recParticles = c12.getDetParticles();   
        int numRec = recParticles.size();     

        bool electronFound = false;  
        float e_vz = -999;
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
                if (!(passesFC && passesVertexCut(pid, vz, e_vz))) {
                    if (i == 0) return;
                    else continue;
                }
                if (det == 1) {
                    float ePCAL  = p->cal(1) ? p->cal(1)->getEnergy() : 0;
                    if (ePCAL <= 0.07) continue;
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
                    e_vz = vz;
                } 
            }
            else if (pid == 2212) {
                bool passesFC = ((det == 1 && FC_->passesDC(p, torus_)) || 
                                 (det == 2 && FC_->passesCVT(p)));
                if (!(passesFC && passesVertexCut(pid, vz, e_vz))) continue;
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
    else if (channel_ == "inclusiveMATCH") {
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

        bool electronFound = false;  
        float e_vz = -999;
        TLorentzVector lv_gen_ePrime;  
        TLorentzVector lv_ePrime; // will store the first electron’s 4-vector

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

            if (gen_pid == 11) {
                // Reject if electron has momentum less than threshold
                if (rec->getP() < 1) continue;

                bool passesFC = ((det == 0 && FC_->passesFT(rec)) || 
                                 (det == 1 && FC_->passesDC(rec, torus_) && FC_->passesECAL(rec)));
                if (!(passesFC && passesVertexCut(gen_pid, vz, e_vz))) continue;

                if (det == 1) {
                    float ePCAL  = rec->cal(1) ? rec->cal(1)->getEnergy() : 0;
                    if (ePCAL <= 0.07) continue;
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
                    e_vz = vz;
                } 
            } 

            else if (gen_pid == 2212) {
                // Reject if proton has momentum less than threshold (Yijie + Bobby)
                if (rec->getP() < 0.3) continue;

                bool passesFC = ((det == 1 && FC_->passesDC(rec, torus_)) || 
                                 (det == 2 && FC_->passesCVT(rec)));
                if (!(passesFC && passesVertexCut(gen_pid, vz, e_vz))) return;
            }

            else if (gen_pid == 22) {
                if (rec->par()->getBeta() < 0.9 || rec->par()->getBeta() > 1.1) continue;
                bool passesFC = (det == 0 && FC_->passesFT(rec)) || (det == 1 && FC_->passesECAL(rec));
                if (!passesFC) return;
            }

            if (electronFound) {
                gen_dis_.fill(lv_gen_ePrime, ebeam_);
                dis_.fill(lv_ePrime, ebeam_);
            }

            ev_.fill(enabledEvBranches_, c12);
            gen_.fill(mcParticles, j);
            
            if (gen_pid == 2212) {
                double p = rec->getP();
                double theta = rec->getTheta();
                double theta_deg = theta * 180.0/M_PI;
                double phi_wrap = rec->getPhi() < 0 ? rec->getPhi() + 2*M_PI : rec->getPhi();

                double p_corr = p + KC_->deltaP(p, theta, det == 1);
                double theta_corr = theta + KC_->deltaTheta(p, theta, det == 1);
                double phi_corr = phi_wrap + + KC_->deltaPhi(p, theta, det == 1);

                if (phi_corr > M_PI) phi_corr -= 2*M_PI;
                if (phi_corr < -M_PI) phi_corr += 2*M_PI;

                rec_.fill(enabledRecBranches_, rec, p_corr, theta_corr, phi_corr);
                // rec_.fill(enabledRecBranches_, rec);
            }
            else rec_.fill(enabledRecBranches_, rec);
            tree_->Fill();
            numFills_++;
        }
    }
    else if (channel_ == "eppi0REC") processEPPI0REC(c12);
    else if (channel_ == "eppi0GEMC") processEPPI0GEMC(c12);
    else if (channel_ == "eppi0MATCH") processEPPI0MATCH(c12);

}

void ProcessManager::processEPPI0REC(clas12::clas12reader& c12) {
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

    double e_vz = -999;
    double leastDeltaM = 999;

    for (const auto& ele : electrons) {

        // Reject if electron has momentum less than threshold
        if (ele->getP() < 1) continue;

        // Reject if electron vertex isn't in target fiducial volume
        if (!passesVertexCut(ele->par()->getPid(), ele->par()->getVz(), e_vz)) continue;

        int detEle = getDetector(ele->par()->getStatus());
        int secEle = ele->getSector();
        // Reject if electron fails fiducial cuts (ONLY if cuts are listed!)
        if (detEle == 0 && !FC_->passesFT(ele)) continue;
        if (detEle == 1 && (!FC_->passesDC(ele, torus_) || !FC_->passesECAL(ele))) continue;

        if (detEle == 1) {
            float ePCAL  = ele->cal(1) ? ele->cal(1)->getEnergy() : 0;
            // Reject if electron fails threshold energy deposit (RGA analysis note)
            if (ePCAL <= 0.07) continue;
            float eECIN  = ele->cal(4) ? ele->cal(4)->getEnergy() : 0;
            float eECOUT = ele->cal(7) ? ele->cal(7)->getEnergy() : 0;
            float sf = (ePCAL + eECIN + eECOUT) / ele->getP();
            if (!SF_->passTriangleCut(ePCAL, eECIN, ele->getP())) continue;
            if (!SF_->passSigmaCut(secEle, sf, ele->getP())) continue;
        }
        // Electron 4-vector
        TLorentzVector lv_e;
        lv_e.SetPxPyPzE(ele->par()->getPx(), ele->par()->getPy(), ele->par()->getPz(),
                        std::sqrt(ELECTRON_MASS * ELECTRON_MASS + ele->getP() * ele->getP()));

        dis_.fill(lv_e, ebeam_);
        e_vz = ele->par()->getVz();

        // Reject if Q2 < 1, W < 2, or y > 0.8, i.e., not in standard DIS region
        if (!channelCheck(dis_.Q2, dis_.W, dis_.y)) continue;

        for (const auto& pro : protons) {

            // Reject if proton has momentum less than threshold (Yijie + Bobby)
            if (pro->getP() < 0.3) continue;

            // Reject if proton vertex isn't in target fiducial volume
            if (!passesVertexCut(pro->par()->getPid(), pro->par()->getVz(), e_vz)) continue;

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

                    // Reject if photons are not detected in same sector (ANDREY):
                    // if (sec1 != sec2) continue;

                    // Reject if photons are detected in same sector as electron (ANDREY):
                    if (sec1 == secEle || sec2 == secEle) continue;

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
                    if (deltaM > 0.15) continue;

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
    
    double p = best_proton->getP();
    double theta = best_proton->getTheta();
    double theta_deg = theta * 180.0/M_PI;
    double phi_wrap = best_proton->getPhi() < 0 ? best_proton->getPhi() + 2*M_PI : best_proton->getPhi();

    double p_corr = p + KC_->deltaP(p, theta, detBestPro == 1);
    double theta_corr = theta + KC_->deltaTheta(p, theta, detBestPro == 1);
    double phi_corr = phi_wrap + KC_->deltaPhi(p, theta, detBestPro == 1);

    if (phi_corr > M_PI) phi_corr -= 2*M_PI;
    if (phi_corr < -M_PI) phi_corr += 2*M_PI;
    
    p_.fill(enabledProBranches_, best_proton, p_corr, theta_corr, phi_corr);
    //p_.fill(enabledProBranches_, best_proton);

    // PHOTON INFO:
    g_.fill(enabledPhoBranches_, best_g1);
    
    // ----- DIS + EPPI0 INFO: ----- //
    
    // Electron 4-vector
    TLorentzVector lv_ePrime;
    lv_ePrime.SetPxPyPzE(best_electron->par()->getPx(),
                         best_electron->par()->getPy(),
                         best_electron->par()->getPz(),
                         std::sqrt(ELECTRON_MASS * ELECTRON_MASS + best_electron->getP() * best_electron->getP()));

    TLorentzVector lv_pPrime;
    lv_pPrime.SetPxPyPzE(p_.p * std::sin(p_.theta) * std::cos(p_.phi),
                         p_.p * std::sin(p_.theta) * std::sin(p_.phi),
                         p_.p * std::cos(p_.theta),
                         std::sqrt(PROTON_MASS * PROTON_MASS + p_.p * p_.p));
    
    dis_.fill(lv_ePrime, ebeam_);

    lv_g1.SetXYZM(best_g1->par()->getPx(), best_g1->par()->getPy(), best_g1->par()->getPz(), 0.0);
    lv_g2.SetXYZM(best_g2->par()->getPx(), best_g2->par()->getPy(), best_g2->par()->getPz(), 0.0);

    eppi0_.fill(lv_ePrime, lv_pPrime, lv_g1, lv_g2, ebeam_);

    // LOOSE GLOBAL EXCLUSIVITY CUTS:
    // if (abs(eppi0_.pT_miss) > 0.2) return;
    // if (eppi0_.t > 2) return;
    if (abs(eppi0_.E_miss) > 2) return;
    if (eppi0_.theta_e_g1 * 180.0/M_PI < 4 || eppi0_.theta_e_g2 * 180.0/M_PI < 4) return; // BOBBY
    if (eppi0_.theta_g1_g2 * 180.0/M_PI < 1) return; // BOBBY
    if (eppi0_.pi0_thetaX * 180.0/M_PI > 4) return; // 2 deg used in CLAS6 analysis
    // if (eppi0_.m2_epX > 1) return;
    // if (eppi0_.m2_epi0X > 3) return;

    numFills_++;
    // FILL TREE, ONCE PER EVENT:
    tree_->Fill();
}

void ProcessManager::processEPPI0GEMC(clas12::clas12reader& c12) { 

    // ===================================================================== //
    // ============= AAO_norad: 1 electron, 1 proton, 2 photons ============ //
    // ============= AAO_rad:   1 electron, 1 proton, 1 pion, 1 photon ===== //

    // =================== GENERATED BLOCK (ALWAYS FILL) =================== //
    TLorentzVector lv_gen_ePrime, lv_gen_pPrime, lv_gen_g1, lv_gen_g2, lv_gen_pi0;
    std::vector<TLorentzVector> lv_genPhotons;
    bool found_pi0 = false;

    auto mcParticles = c12.mcparts();
    int numGen = mcParticles->getRows(); 

    // --------------------------
    // Identify radiative events
    // --------------------------
    bool is_rad = false;
    for (int j = 0; j < numGen; j++) {
        if (mcParticles->getPid(j) == 111) {   // π0 in generator → radiative sample
            is_rad = true;
            break;
        }
    }

    for (int j = 0; j < numGen; j++) {
            int gen_pid = mcParticles->getPid(j);
            TVector3 v_gen_p(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j));    

            if (gen_pid == 11) {
                lv_gen_ePrime.SetPxPyPzE(v_gen_p.Px(), v_gen_p.Py(), v_gen_p.Pz(), 
                                         std::sqrt(ELECTRON_MASS * ELECTRON_MASS + v_gen_p.Mag() * v_gen_p.Mag()));
            }
            else if (gen_pid == 2212) {
                lv_gen_pPrime.SetPxPyPzE(v_gen_p.Px(), v_gen_p.Py(), v_gen_p.Pz(), 
                                         std::sqrt(PROTON_MASS * PROTON_MASS + v_gen_p.Mag() * v_gen_p.Mag()));
            }
            else if (gen_pid == 22) {
                lv_genPhotons.emplace_back(v_gen_p.Px(), v_gen_p.Py(), v_gen_p.Pz(), v_gen_p.Mag());
            }
            else if (gen_pid == 111) {
                lv_gen_pi0.SetPxPyPzE(v_gen_p.Px(), v_gen_p.Py(), v_gen_p.Pz(), 
                                         std::sqrt(PI0_MASS * PI0_MASS + v_gen_p.Mag() * v_gen_p.Mag()));
            }
    }

    if (!is_rad) {
        // NONRADIATIVE: expect 2 generator photons
        if (lv_genPhotons.size() < 2) {
            // you may want to print a warning here
            return;
        }
        lv_gen_g1 = lv_genPhotons[0];
        lv_gen_g2 = lv_genPhotons[1];
    }
    else {
        // RADIATIVE: expect π0 + radiative γ
        // Generator does not give π0 → γγ decay products (GEANT does)
        // The radiative photon is the single generator γ
        if (lv_genPhotons.size() != 1) {
            // RadGen should produce exactly one photon at generator level
            // If not: skip or debug warning.
            return;
        } 
        lv_gen_g1 = lv_genPhotons[0]; // assign radiative photon to g1
    }
    
    gen_dis_.flush();
    gen_dis_.fill(lv_gen_ePrime, ebeam_);

    gen_eppi0_.flush();

    // Choose objects to pass in based on mode
    TLorentzVector lv_obj1 = is_rad ? lv_gen_pi0 : lv_gen_g1;
    TLorentzVector lv_obj2 = is_rad ? lv_gen_g1 : lv_gen_g2;

    gen_eppi0_.fill(lv_gen_ePrime, lv_gen_pPrime, lv_obj1, lv_obj2, ebeam_, is_rad);

    // Reject if Q2 < 1, W < 2, or y > 0.8, i.e., not in standard DIS region  Need to think about removing this for migrations
    if (!channelCheck(gen_dis_.Q2, gen_dis_.W, gen_dis_.y)) return;

    // =================== RECONSTRUCTED BLOCK =================== //

    bool passREC = false;

    auto electrons = c12.getByID(11);   
    auto protons   = c12.getByID(2212); 
    auto photons   = c12.getByID(22); 

    if (photons.size() >= 2 || electrons.size() == 1 || protons.size() == 1) {

        clas12::region_particle* best_electron = nullptr;
        clas12::region_particle* best_proton   = nullptr;
        clas12::region_particle* best_g1       = nullptr;
        clas12::region_particle* best_g2       = nullptr;

        TLorentzVector lv_g1, lv_g2;
        int detPi0;
        double e_vz = -999;
        double leastDeltaM = 999;

        for (const auto& ele : electrons) {

            if (ele->getP() < 1) continue;
            if (!passesVertexCut(ele->par()->getPid(), ele->par()->getVz(), e_vz)) continue;

            int detEle = getDetector(ele->par()->getStatus());
            int secEle = ele->getSector();
        
            if (detEle == 0 && !FC_->passesFT(ele)) continue;
            if (detEle == 1 && (!FC_->passesDC(ele, torus_) || !FC_->passesECAL(ele))) continue;

            if (detEle == 1) {
                float ePCAL  = ele->cal(1) ? ele->cal(1)->getEnergy() : 0;
                if (ePCAL <= 0.07) continue;
                float eECIN  = ele->cal(4) ? ele->cal(4)->getEnergy() : 0;
                float eECOUT = ele->cal(7) ? ele->cal(7)->getEnergy() : 0;
                float sf = (ePCAL + eECIN + eECOUT) / ele->getP();
                if (!SF_->passTriangleCut(ePCAL, eECIN, ele->getP())) continue;
                if (!SF_->passSigmaCut(secEle, sf, ele->getP())) continue;
            }
            
            TLorentzVector lv_e;
            lv_e.SetPxPyPzE(ele->par()->getPx(), ele->par()->getPy(), ele->par()->getPz(),
                            std::sqrt(ELECTRON_MASS * ELECTRON_MASS + ele->getP() * ele->getP()));

            dis_.fill(lv_e, ebeam_);
            e_vz = ele->par()->getVz();

            if (!channelCheck(dis_.Q2, dis_.W, dis_.y)) continue;

            for (const auto& pro : protons) {

                if (pro->getP() < 0.3) continue;
                if (!passesVertexCut(pro->par()->getPid(), pro->par()->getVz(), e_vz)) continue;

                int detPro = getDetector(pro->par()->getStatus());
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

                        if (det1 == 2 || det2 == 2) continue;
                        //if (sec1 != sec2) continue;

                        if (sec1 == secEle || sec2 == secEle) continue;

                        if (g1->getP() < 0.4  || g2->getP() < 0.4) continue;

                        if (g1->par()->getBeta() < 0.9 || g1->par()->getBeta() > 1.1) continue;
                        if (g2->par()->getBeta() < 0.9 || g2->par()->getBeta() > 1.1) continue;

                        if (getPhotonCalE(g1) < 0.15 || getPhotonCalE(g2) < 0.15) continue;

                        if ((det1 == 0 && !FC_->passesFT(g1)) || (det2 == 0 && !FC_->passesFT(g2))) continue;
                        if ((det1 == 1 && !FC_->passesECAL(g1)) || (det2 == 1 && !FC_->passesECAL(g2))) continue;

                        lv_g1.SetXYZM(g1->par()->getPx(), g1->par()->getPy(), g1->par()->getPz(), 0.0);
                        lv_g2.SetXYZM(g2->par()->getPx(), g2->par()->getPy(), g2->par()->getPz(), 0.0);
                        TLorentzVector lv_candidatePi0 = lv_g1 + lv_g2;

                        double deltaM = std::abs(lv_candidatePi0.M() - PI0_MASS);
                        if (deltaM > 0.15) continue;

                        if (deltaM < leastDeltaM) {
                            detPi0         = det1;
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

        if (best_electron) {
            ev_.fill(enabledEvBranches_, c12);

            e_.fill(enabledEleBranches_, best_electron);

            int detBestPro = getDetector(best_proton->par()->getStatus());
            double p = best_proton->getP();
            double theta = best_proton->getTheta();
            double theta_deg = theta * 180.0/M_PI;
            double phi_wrap = best_proton->getPhi() < 0 ? best_proton->getPhi() + 2*M_PI : best_proton->getPhi();
            double p_corr = p + KC_->deltaP(p, theta, detBestPro == 1);
            double theta_corr = theta + KC_->deltaTheta(p, theta, detBestPro == 1);
            double phi_corr = phi_wrap + KC_->deltaPhi(p, theta, detBestPro == 1);
            if (phi_corr > M_PI) phi_corr -= 2*M_PI;
            if (phi_corr < -M_PI) phi_corr += 2*M_PI;
            
            p_.fill(enabledProBranches_, best_proton, p_corr, theta_corr, phi_corr);

            // PHOTON INFO:
            g_.fill(enabledPhoBranches_, best_g1);
            
            // ----- DIS + EPPI0 INFO: ----- //
            
            // Electron 4-vector
            TLorentzVector lv_ePrime;
            lv_ePrime.SetPxPyPzE(best_electron->par()->getPx(),
                                best_electron->par()->getPy(),
                                best_electron->par()->getPz(),
                                std::sqrt(ELECTRON_MASS * ELECTRON_MASS + best_electron->getP() * best_electron->getP()));
            TLorentzVector lv_pPrime;
            lv_pPrime.SetPxPyPzE(p_.p * std::sin(p_.theta) * std::cos(p_.phi),
                                p_.p * std::sin(p_.theta) * std::sin(p_.phi),
                                p_.p * std::cos(p_.theta),
                                std::sqrt(PROTON_MASS * PROTON_MASS + p_.p * p_.p));
            dis_.fill(lv_ePrime, ebeam_);

            lv_g1.SetXYZM(best_g1->par()->getPx(), best_g1->par()->getPy(), best_g1->par()->getPz(), 0.0);
            lv_g2.SetXYZM(best_g2->par()->getPx(), best_g2->par()->getPy(), best_g2->par()->getPz(), 0.0);

            eppi0_.fill(lv_ePrime, lv_pPrime, lv_g1, lv_g2, ebeam_);

            do {
                // LOOSE GLOBAL EXCLUSIVITY CUTS:
                // if (eppi0_.pT_miss > 0.2) return;
                // if (eppi0_.t > 2) return;
                if (abs(eppi0_.E_miss) > 2) return;
                if (eppi0_.theta_e_g1 * 180.0/M_PI < 4 || eppi0_.theta_e_g2 * 180.0/M_PI < 4) return; // BOBBY
                if (eppi0_.theta_g1_g2 * 180.0/M_PI < 1) return; // BOBBY
                if (eppi0_.pi0_thetaX * 180.0/M_PI > 4) return; // used in CLAS6 analysis
                // if (eppi0_.m2_epX > 1) return;
                // if (eppi0_.m2_epi0X > 3) return;

                passREC = true;

            } while(0);
        }
    }
    if (!passREC) {
        ev_.flush();
        e_.flush();
        p_.flush();
        g_.flush();
        dis_.flush();
        eppi0_.flush();
    }
    numFills_++;
    tree_->Fill();
}

void ProcessManager::processEPPI0MATCH(clas12::clas12reader& c12) {
    auto electrons = c12.getByID(11);   // PID 11 = electron
    auto protons   = c12.getByID(2212); // PID 2212 = proton
    auto photons   = c12.getByID(22);   // PID 22 = photon

    // Reject if not EPPI0 final state, according to eventbuilder:
    if (photons.size() < 2 || electrons.size() != 1 || protons.size() != 1) return;

    TLorentzVector lv_gen_ePrime, lv_gen_pPrime, lv_gen_g1, lv_gen_g2;
    TLorentzVector lv_ePrime, lv_pPrime, lv_g1, lv_g2;
    std::vector<std::pair<int,int>> photonMatches;  // (gen indx, rec indx)

    auto mcParticles = c12.mcparts();
    auto recParticles = c12.getDetParticles();

    int numGen = mcParticles->getRows(); 
    int numRec = recParticles.size();

    float e_vz = -999;

    for (int j = 0; j < numGen; j++) {
            int gen_pid = mcParticles->getPid(j);
            TVector3 v_gen_p(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j));

            double bestDeltaTheta = 999.0;
            int bestMatchIndex = -1;

            for (int i = 0; i < numRec; i++) {
                auto& rec = recParticles[i];
                if (rec->par()->getPid() != gen_pid) continue;

                TVector3 v_rec_p(rec->par()->getPx(), rec->par()->getPy(), rec->par()->getPz());
                double deltaTheta = v_gen_p.Angle(v_rec_p) * 180/M_PI; 

                if (deltaTheta < bestDeltaTheta) {
                    bestDeltaTheta = deltaTheta;
                    bestMatchIndex = i;
                }
            }

            // Threshold for match quality:
            if (bestMatchIndex < 0 || bestDeltaTheta > 3) {
                if (gen_pid == 11 || gen_pid == 2212) return;
                continue;
            }

            auto& rec  = recParticles[bestMatchIndex];
            int det    = getDetector(rec->par()->getStatus());
            int sec    = rec->getSector();
            float vz   = rec->par()->getVz();      

            if (gen_pid == 11) {
                // Reject if electron has momentum less than threshold
                if (rec->getP() < 1) return;

                bool passesFC = ((det == 0 && FC_->passesFT(rec)) || 
                                 (det == 1 && FC_->passesDC(rec, torus_) && FC_->passesECAL(rec)));
                if (!(passesFC && passesVertexCut(gen_pid, vz, e_vz))) return;

                if (det == 1) {
                    float ePCAL  = rec->cal(1) ? rec->cal(1)->getEnergy() : 0;
                    if (ePCAL <= 0.07) continue;
                    float eECIN  = rec->cal(4) ? rec->cal(4)->getEnergy() : 0;
                    float eECOUT = rec->cal(7) ? rec->cal(7)->getEnergy() : 0;
                    float sf = (ePCAL + eECIN + eECOUT) / rec->getP();
                    if (!SF_->passTriangleCut(ePCAL, eECIN, rec->getP())) return;
                    if (!SF_->passSigmaCut(sec, sf, rec->getP())) return;
                }

                lv_gen_ePrime.SetPxPyPzE(v_gen_p.Px(), 
                                         v_gen_p.Py(), 
                                         v_gen_p.Pz(), 
                                         std::sqrt(ELECTRON_MASS * ELECTRON_MASS + v_gen_p.Mag() * v_gen_p.Mag()));

                lv_ePrime.SetPxPyPzE(rec->par()->getPx(), 
                                     rec->par()->getPy(), 
                                     rec->par()->getPz(), 
                                     std::sqrt(ELECTRON_MASS * ELECTRON_MASS + rec->getP() * rec->getP()));

                dis_.flush();
                dis_.fill(lv_ePrime, ebeam_);
                // Reject if Q2 < 1, W < 2, or y > 0.8, i.e., not in standard DIS region
                if (!channelCheck(dis_.Q2, dis_.W, dis_.y)) return;

                gen_dis_.flush();
                gen_dis_.fill(lv_gen_ePrime, ebeam_);

                e_.flush();
                e_.fill(enabledEleBranches_, rec);
                e_vz = vz;
            }

            else if (gen_pid == 2212) {
                // Reject if proton has momentum less than threshold (Yijie + Bobby)
                if (rec->getP() < 0.3) return;
                bool passesFC = ((det == 1 && FC_->passesDC(rec, torus_)) || 
                                 (det == 2 && FC_->passesCVT(rec)));
                if (!(passesFC && passesVertexCut(gen_pid, vz, e_vz))) return;

                double p = rec->getP();
                double theta = rec->getTheta();
                double theta_deg = theta * 180.0/M_PI;
                double phi_wrap = rec->getPhi() < 0 ? rec->getPhi() + 2*M_PI : rec->getPhi();

                double p_corr = p + KC_->deltaP(p, theta, det == 1);
                double theta_corr = theta + KC_->deltaTheta(p, theta, det == 1);
                double phi_corr = phi_wrap + + KC_->deltaPhi(p, theta, det == 1);

                if (phi_corr > M_PI) phi_corr -= 2*M_PI;
                if (phi_corr < -M_PI) phi_corr += 2*M_PI;

                p_.flush();
                p_.fill(enabledProBranches_, rec, p_corr, theta_corr, phi_corr);

                lv_gen_pPrime.SetPxPyPzE(v_gen_p.Px(), 
                                         v_gen_p.Py(), 
                                         v_gen_p.Pz(), 
                                         std::sqrt(PROTON_MASS * PROTON_MASS + v_gen_p.Mag() * v_gen_p.Mag()));

                lv_pPrime.SetPxPyPzE(p_.p * std::sin(p_.theta) * std::cos(p_.phi),
                                     p_.p * std::sin(p_.theta) * std::sin(p_.phi),
                                     p_.p * std::cos(p_.theta),
                                     std::sqrt(PROTON_MASS * PROTON_MASS + p_.p * p_.p));
            }

            else if (gen_pid == 22) {
                if (det == 2) continue;
                if (rec->getP() < 0.4) continue;
                if (rec->par()->getBeta() < 0.9 || rec->par()->getBeta() > 1.1) continue;
                if (getPhotonCalE(rec) < 0.15) continue;
                bool passesFC = (det == 0 && FC_->passesFT(rec)) || (det == 1 && FC_->passesECAL(rec));
                if (!passesFC) continue;
                photonMatches.emplace_back(j, bestMatchIndex); // save GEN and REC indices;

                g_.flush();
                g_.fill(enabledPhoBranches_, rec);
            }
    }

    // now test photon pairings
    if (photonMatches.size() < 2) return;
    
    bool foundPi0 = false;
    double leastDeltaM = 1e9;

    for (size_t i = 0; i < photonMatches.size(); i++) {
        for (size_t j = i+1; j < photonMatches.size(); j++) {

            auto& g1 = recParticles[photonMatches[i].second];
            auto& g2 = recParticles[photonMatches[j].second];

            int det1 = getDetector(g1->par()->getStatus());
            int det2 = getDetector(g2->par()->getStatus());

            int sec1 = g1->getSector();
            int sec2 = g2->getSector();

            // Reject if photons are not detected in same sector
            if (sec1 != sec2) continue;

            if (sec1 == e_.sector || sec2 == e_.sector) continue;

            TLorentzVector lv_tmp_g1, lv_tmp_g2, lv_candidatePi0;
            lv_tmp_g1.SetXYZM(g1->par()->getPx(), g1->par()->getPy(), g1->par()->getPz(), 0.0);
            lv_tmp_g2.SetXYZM(g2->par()->getPx(), g2->par()->getPy(), g2->par()->getPz(), 0.0);
            lv_candidatePi0 = lv_tmp_g1 + lv_tmp_g2;

            double deltaM = std::abs(lv_candidatePi0.M() - PI0_MASS);
            if (deltaM > 0.1) continue;

            if (deltaM < leastDeltaM) {
                foundPi0 = true;
                leastDeltaM = deltaM;
                int gen_g1 = photonMatches[i].first;
                int gen_g2 = photonMatches[j].first;
                lv_gen_g1.SetXYZM(mcParticles->getPx(gen_g1), mcParticles->getPy(gen_g1), mcParticles->getPz(gen_g1), 0);
                lv_gen_g2.SetXYZM(mcParticles->getPx(gen_g2), mcParticles->getPy(gen_g2), mcParticles->getPz(gen_g2), 0);
                lv_g1       = lv_tmp_g1;
                lv_g2       = lv_tmp_g2;
            }
        }
    }

    if (!foundPi0) return;
    
    gen_eppi0_.flush();
    eppi0_.flush();

    gen_eppi0_.fill(lv_gen_ePrime, lv_gen_pPrime, lv_gen_g1, lv_gen_g2, ebeam_);
    eppi0_.fill(lv_ePrime, lv_pPrime, lv_g1, lv_g2, ebeam_);

    // LOOSE GLOBAL EXCLUSIVITY CUTS:
    if (eppi0_.m2_miss > 1) return;
    if (eppi0_.pT_miss > 0.2) return;
    if (eppi0_.t > 2) return;
    if (eppi0_.E_miss > 1.5) return;
    if (eppi0_.pz_miss > 1) return;
    if (eppi0_.theta_e_g1 * 180.0/M_PI < 4 || eppi0_.theta_e_g2 * 180.0/M_PI < 4) return; // BOBBY
    if (eppi0_.theta_g1_g2 * 180.0/M_PI < 1) return; // BOBBY
    if (eppi0_.pi0_thetaX * 180.0/M_PI > 4) return; // used in CLAS6 analysis
    if (eppi0_.m2_epX > 1) return;
    if (eppi0_.m2_epi0X > 3) return;

    numFills_++;

    // FILL TREE, ONCE PER EVENT:
    tree_->Fill();

}
