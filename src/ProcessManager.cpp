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

    ebeam_       = config["ebeam"];
    torus_       = config["torus"];
    channel_     = config["channel"];
    outPrefix_   = config.value("outPrefix", "");

    FC_ = new FiducialCuts();

    for (const auto& cut : config.value("fiducialCuts", std::vector<std::string>{})) { FC_->addCut(cut); }

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
    tree_->Branch("DIS", "DISVars", &dis_);

    if (channel_ == "aaoGenOnly") {
        tree_->Branch("gen", "GenVars", &gen_);
    }

    else if (channel_ == "inclusiveRec") {
        tree_->Branch("rec", "RecVars", &rec_);
    }

    else if (channel_ == "aaoGenMatch") {
        tree_->Branch("GEN_DIS",  "DISVars",     &gen_dis_);
        tree_->Branch("gen",      "GenVars",     &gen_);
        tree_->Branch("rec",      "RecVars",     &rec_);
    }

    else if (channel_ == "eppi0") {
        tree_->Branch("e",     "RecVars",     &e_);
        tree_->Branch("p",     "RecVars",     &p_);
        tree_->Branch("g",     "RecVars",     &g_);
        tree_->Branch("eppi0", "EPPI0Vars",   &eppi0_);
    }

    if (config.contains("branches") && !config["branches"].is_null()) {
        auto& b = config["branches"];
        for (const auto& var : b.value("event", std::vector<std::string>{}))
            enabledEvBranches_.insert(var);
        for (const auto& var : b.value("rec", std::vector<std::string>{}))
            enabledRecBranches_.insert(var);
        for (const auto& var : b.value("electron", std::vector<std::string>{}))
            enabledEleBranches_.insert(var);
        for (const auto& var : b.value("proton", std::vector<std::string>{}))
            enabledProBranches_.insert(var);
        for (const auto& var : b.value("photon", std::vector<std::string>{}))
            enabledPhoBranches_.insert(var);
    } 
}

//////////////////////////////////////////////////////////
// Little helper functions called by bigger functions: ///
//////////////////////////////////////////////////////////

std::string ProcessManager::currentTimestamp() const {
    auto now = std::chrono::system_clock::now();
    std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&tt), "%H%M_%m%d%y");
    return ss.str();
}

// based on particle status, returns int representation of FT (0), FD (1), or CD (2)
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

//////////////////////////////////////////////////////////
// Bigger functions: //
//////////////////////////////////////////////////////////

void ProcessManager::finalize() {
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
    eventsProcessed_++;
    if (channel_ == "aaoGenOnly") {

        // Generated particles
        auto mcParticles = c12.mcparts();
        int numGen = mcParticles->getRows();       

        bool electronFilled = false;  
        TLorentzVector lv_ePrime; // will store the first electron’s 4-vector

        for (int j = 0; j < numGen; j++) {

            ev_.fill(enabledEvBranches_, c12);
            gen_.fill(mcParticles, j);

            int pid = mcParticles->getPid(j);
            TVector3 v_p(mcParticles->getPx(j), mcParticles->getPy(j), mcParticles->getPz(j));

            if (pid == 11 && !electronFilled) {
                // Calculate and store only for the first electron
                lv_ePrime.SetVectM(v_p, ELECTRON_MASS);
                dis_.fill(lv_ePrime, ebeam_);
                electronFilled = true;
            } 

            else if (!electronFilled) {
                dis_.flush();
            }

            else dis_.fill(lv_ePrime, ebeam_);

            tree_->Fill();
            numFills_++;
        }

    }
    else if (channel_ == "inclusiveRec") {
        
        auto recParticles = c12.getDetParticles();   
        int numRec = recParticles.size();     

        // Initialize DIS vars
        dis_.flush();
        bool electronFilled = false;  
        TLorentzVector lv_ePrime; // will store the first electron’s 4-vector       

        for (int i = 0; i < numRec; i++) {
            auto& p = recParticles[i];

            int pid  = p->getPid();
            int det  = getDetector(p->par()->getStatus());
            float vz = p->par()->getVz();

            if (pid == 11) {
                bool passesFC = ((det == 0 && FC_->passesFT(p)) || 
                                 (det == 1 && FC_->passesDC(p, torus_) && FC_->passesECAL(p)));
                if (!(passesFC && passesVertexCut(vz))) return;

                if (!electronFilled) {
                    lv_ePrime.SetPxPyPzE(p->par()->getPx(), 
                                         p->par()->getPy(), 
                                         p->par()->getPz(), 
                                         ELECTRON_MASS * ELECTRON_MASS + p->getP() * p->getP());
                    dis_.fill(lv_ePrime, ebeam_);
                    electronFilled = true;
                } 
            }

            else if (pid == 2212) {
                bool passesFC = ((det == 1 && FC_->passesDC(p, torus_)) || 
                                 (det == 2 && FC_->passesCVT(p)));
                if (!(passesFC && passesVertexCut(vz))) return;
            }
            else if (pid == 22) {
                bool passesFC = (det == 0 && FC_->passesFT(p)) || (det == 1 && FC_->passesECAL(p));
                if (!passesFC) return;
            }

            // If electron already found, repeat DIS fill for all subsequent particles
            if (electronFilled) dis_.fill(lv_ePrime, ebeam_);

            ev_.fill(enabledEvBranches_, c12);
            rec_.fill(enabledRecBranches_, p);
            tree_->Fill();
            numFills_++;
        }
    }
    else if (channel_ == "aaoGenMatch") {
        
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
            float vz   = rec->par()->getVz();

            // Initialize DIS vars
            dis_.flush();
            gen_dis_.flush();
            bool electronFilled = false;  
            TLorentzVector lv_gen_ePrime;  
            TLorentzVector lv_ePrime; // will store the first electron’s 4-vector       

            if (gen_pid == 11) {
                bool passesFC = ((det == 0 && FC_->passesFT(rec)) || 
                                 (det == 1 && FC_->passesDC(rec, torus_) && FC_->passesECAL(rec)));
                if (!(passesFC && passesVertexCut(vz))) return;

                if (!electronFilled) {
                    lv_gen_ePrime.SetPxPyPzE(v_gen_p.Px(), 
                                         v_gen_p.Py(), 
                                         v_gen_p.Pz(), 
                                         ELECTRON_MASS * ELECTRON_MASS + v_gen_p.Mag() * v_gen_p.Mag());

                    gen_dis_.fill(lv_gen_ePrime, ebeam_);

                    lv_ePrime.SetPxPyPzE(rec->par()->getPx(), 
                                         rec->par()->getPy(), 
                                         rec->par()->getPz(), 
                                         ELECTRON_MASS * ELECTRON_MASS + rec->getP() * rec->getP());
                    dis_.fill(lv_ePrime, ebeam_);
                    electronFilled = true;
                } 
            } 

            else if (gen_pid == 2212) {
                bool passesFC = ((det == 1 && FC_->passesDC(rec, torus_)) || 
                                 (det == 2 && FC_->passesCVT(rec)));
                if (!(passesFC && passesVertexCut(vz))) return;
            }

            else if (gen_pid == 22) {
                bool passesFC = (det == 0 && FC_->passesFT(rec)) || (det == 1 && FC_->passesECAL(rec));
                if (!passesFC) return;
            }

            // If electron already found, repeat DIS fill for all subsequent particles
            if (electronFilled) {
                gen_dis_.fill(lv_gen_ePrime, ebeam_);
                dis_.fill(lv_ePrime, ebeam_);
            }

            ev_.fill(enabledEvBranches_, c12);
            gen_.fill(mcParticles, j);
            rec_.fill(enabledRecBranches_, rec);
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

    clas12::region_particle* best_electron = nullptr;
    clas12::region_particle* best_proton   = nullptr;
    clas12::region_particle* best_g1       = nullptr;
    clas12::region_particle* best_g2       = nullptr;
    
    TLorentzVector lv_bestPi0;
    int detPi0;

    double  leastDeltaM = 999;
    double  m_bestPi0 = 999;

    for (const auto& ele : electrons) {

        // Reject if electron vertex isn't in target fiducial volume
        if (!passesVertexCut(ele->par()->getVz())) continue;

        int detEle = getDetector(ele->par()->getStatus());
        // Reject if electron fails fiducial cuts (ONLY if cuts are listed!)
        if (detEle == 0 && !FC_->passesFT(ele)) continue;
        if (detEle == 1 && (!FC_->passesDC(ele, torus_) || !FC_->passesECAL(ele))) continue;

        for (const auto& pro : protons) {

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

                    // Reject if either photon is detected in CD
                    if (det1 == 2 || det2 == 2) continue;

                    int sec1 = g1->getSector();
                    int sec2 = g2->getSector();

                    // Reject if photons are not detected in same sector
                    if (sec1 != sec2) continue;

                    // Reject if either photon is detected in FT and fails FT fiducial cuts (fails ONLY if FTstandardCut is listed!)
                    if ((det1 == 0 && !FC_->passesFT(g1)) || (det2 == 0 && !FC_->passesFT(g2))) continue;

                    //Reject if either photon is detected in FD and fails ECAL fiducial cuts (fails ONLY if ECAL cuts are listed!)
                    if ((det1 == 1 && !FC_->passesECAL(g1)) || (det2 == 1 && !FC_->passesECAL(g2))) continue;


                    TLorentzVector lv_g1, lv_g2;
                    lv_g1.SetXYZM(g1->par()->getPx(), g1->par()->getPy(), g1->par()->getPz(), 0.0);
                    lv_g2.SetXYZM(g2->par()->getPx(), g2->par()->getPy(), g2->par()->getPz(), 0.0);
                    TLorentzVector lv_candidatePi0 = lv_g1 + lv_g2;

                    double deltaM = std::abs(lv_candidatePi0.M() - PI0_MASS);
                    if (deltaM > 0.025) continue;

                    if (deltaM < leastDeltaM) {
                        m_bestPi0      = lv_candidatePi0.M();
                        detPi0         = det1; // Assign to pi0 the detector of the photon
                        leastDeltaM    = deltaM;
                        lv_bestPi0     = lv_candidatePi0;
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

    ev_.flush();
    e_.flush();
    p_.flush();
    g_.flush();
    dis_.flush();
    eppi0_.flush();

    ev_.fill(enabledEvBranches_, c12);

    // ELECTRON INFO:
    e_.fill(enabledEleBranches_, best_electron);

    // PROTON INFO:
    p_.fill(enabledProBranches_, best_proton);

    // PHOTON INFO:
    g_.fill(enabledPhoBranches_, best_g1);
    
    // DIS + EPPI0 INFO:
    
    // Electron 4-vector
    TLorentzVector lv_ePrime;
    lv_ePrime.SetPxPyPzE(e_.px,e_.py, e_.pz, std::sqrt(ELECTRON_MASS * ELECTRON_MASS + e_.p * e_.p));

    TLorentzVector lv_pPrime;
    lv_pPrime.SetPxPyPzE(p_.px, p_.py, p_.pz, std::sqrt(PROTON_MASS * PROTON_MASS + p_.p * p_.p));
    
    dis_.fill(lv_ePrime, ebeam_);

    eppi0_.fill(lv_ePrime, lv_pPrime, lv_bestPi0, ebeam_);

    // Reject if Q2 < 1, W < 2, or y > 0.8, i.e., not in standard DIS region
    if (!channelCheck(dis_.Q2, dis_.W, dis_.y)) return;

    numFills_++;

    // FILL TREE, ONCE PER EVENT:
    tree_->Fill();
}
