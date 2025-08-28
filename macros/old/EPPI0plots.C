// To run from command line, type << root -l -b -q 'rootMacros/EPPI0plots.C("output/test.root")' >>

// Check if all required branches exist in the tree.
// Returns true if all found, else prints missing and returns false.
bool checkBranches(TTree* tree, const std::vector<std::string>& branches) {
    bool allFound = true;
    std::vector<std::string> missingBranches;
    for (const auto& b : branches) {
        if (!tree->GetBranch(b.c_str())) {
            missingBranches.push_back(b);
            allFound = false;
        }
    }
    if (!allFound) {
        std::cerr << "\nERROR: Missing required branches in the TTree 'Events':\n";
        for (const auto& mb : missingBranches) {
            std::cerr << "  - " << mb << '\n';
        }
        std::cerr << "\nPlease ensure the input ROOT file and TTree contain all these branches.\n"
                     "Required branches:\n";
        for (const auto& b : branches) {
            std::cerr << "  " << b << '\n';
        }
    }
    return allFound;
}

void EPPI0plots(const char* inputFilePath = "input.root", const char* outFilePath = "EPPI0plots.root") {
    // ─── Open Input File and Tree ─────────────────────────────

    gSystem->Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so");


    // List of required branches for this macro to run correctly.
    const std::vector<std::string> requiredBranches = {
        "event.helicity",
        "dis.Q2", "dis.nu", "dis.Xb", "dis.y", "dis.W",
        "e.p", "e.theta", "e.phi", "e.vz", "e.chi2pid", "e.det", "e.sector",
        "e.E_PCAL", "e.E_ECIN", "e.E_ECOUT",
        "e.xDC1", "e.yDC1", "e.xDC2", "e.yDC2", "e.xDC3", "e.yDC3",
        "p.p", "p.theta", "p.phi", "p.vz", "p.chi2pid", "p.det", "p.sector",
        "p.edge_cvt1", "p.edge_cvt3", "p.edge_cvt5", "p.edge_cvt7", "p.edge_cvt12",
        "p.theta_cvt", "p.phi_cvt",
        "p.xDC1", "p.yDC1", "p.xDC2", "p.yDC2", "p.xDC3", "p.yDC3",
        "g.p", "g.theta", "g.phi", "g.chi2pid", "g.det", "g.sector",
        "g.E_PCAL", "g.E_ECIN", "g.E_ECOUT",
        "e.vPCAL", "e.wPCAL", "e.vECIN", "e.wECIN", "e.vECOUT", "e.wECOUT", "e.uECOUT",
        "g.vPCAL", "g.wPCAL", "g.vECIN", "g.wECIN", "g.vECOUT", "g.wECOUT", "g.uECOUT"
    };

    TFile* file = TFile::Open(inputFilePath);
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open file: " << inputFilePath << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "Could not find tree 'Events' in file." << std::endl;
        return;
    }

    if (!checkBranches(tree, requiredBranches)) {
        // Missing required branches, abort.
        std::cout << "Proceeding with remaining histograms." << std::endl;
    }

    // ─── Create Histograms ───────────────────────────────────
    TH1I* h_helicity = new TH1I("h_helicity", "Helicity; Counts", 10, -5, 5);
    //TH1I* h_numPhotons = new TH1I("h_numPhotons", "# of Photons; Counts", 6, 0, 6);
    TH1D* h_mgg = new TH1D("h_mgg", "M_{#gamma#gamma}; M_{(#gamma#gamma)} [GeV]; Counts", 200, 0.098, 0.17);
    TH1D* h_Q2 = new TH1D("h_Q2", "Q^{2} Distribution; Q^{2} [GeV^{2}]; Counts", 250, 0, 10);
    TH1D* h_W = new TH1D("h_W", "W Distribution; W [GeV]; Counts", 100, 2, 3.5);
    TH1D* h_t = new TH1D("h_t", "-t Distribution; -t [GeV^{2}]; Counts", 100, 0, 4);
    TH1D* h_phiT = new TH1D("h_phiT", "Trento #Phi Distribution; #Phi_{T} [deg]; Counts", 180, 0, 360);
    TH1D* h_m2_epX = new TH1D("h_m2_epX", "M_{epX}^{2}; M_{epX}^{2} [GeV^{2}]; Counts", 200, -1, 4);
    TH1D* h_m2_epi0X = new TH1D("h_m2_epi0X", "M_{e#piX}^{2}; M_{e#piX}^{2} [GeV^{2}]; Counts", 200, -1, 4);
    TH1D* h_px_miss = new TH1D("h_px_miss", "Missing Momentum #Delta P_{x}; #Delta P_{x} [GeV]; Counts", 200, -0.8, 0.8);
    TH1D* h_py_miss = new TH1D("h_py_miss", "Missing Momentum #Delta P_{y}; #Delta P_{y} [GeV]; Counts", 200, -0.8, 0.8);
    TH1D* h_pz_miss = new TH1D("h_pz_miss", "Missing Momentum #Delta P_{z}; #Delta P_{z} [GeV]; Counts", 200, -0.3, 0.6);
    TH1D* h_E_miss = new TH1D("h_E_miss", "E_{miss}; E_{miss} [GeV]; Counts", 200, -0.8, 1.2);
    TH1D* h_deltaPhi = new TH1D("h_deltaPhi", "#Delta #Phi #equiv #Phi _{#pi} - #Phi _{epX}; #Delta #Phi [Deg]; Counts", 200, -10, 10);

    TH2D* h_Q2_vs_Xb = new TH2D("h_Q2_vs_Xb", "Q^{2} vs x_{B}; x_{B}; Q^{2} [GeV^{2}]", 
                        50, 0.0, 0.7, 50, 0.5, 7);
    TH2D* h_t_vs_phiT = new TH2D("h_t_vs_phiT", "-t vs. #Phi_{T}; #Phi_{T} [deg]; -t [Gev^{2}]", 
                        180, 0, 360, 50, 0, 2);
    TH2D* h_e_theta_vs_p = new TH2D("h_e_theta_vs_p", "Electron #theta vs p; Electron Momentum p [GeV]; #theta [deg]", 
                        50, 1, 4.5, 80, 0, 40);
    TH2D* h_e_phi_vs_p = new TH2D("h_e_phi_vs_p", "Electron #Phi vs p; Electron Momentum p [GeV]; #Phi [deg]", 
                        50, 1, 4.5, 180, -180, 180);
    TH2D* h_p_theta_vs_p = new TH2D("h_p_theta_vs_p", "Proton #theta vs p; Proton Momentum p [GeV]; #theta [deg]", 
                        50, 0, 4.5, 100, 0, 90);
    TH2D* h_p_phi_vs_p = new TH2D("h_p_phi_vs_p", "Proton #Phi vs p; Proton Momentum p [GeV]; #Phi [deg]", 
                        50, 0, 4.5, 100, -180, 180);
    TH2D* h_g_theta_vs_p = new TH2D("h_g_theta_vs_p", "Photon #theta vs p; Photon Momentum p [GeV]; #theta [deg]", 
                        50, 0, 5, 100, 0, 40);
    TH2D* h_g_phi_vs_p = new TH2D("h_g_phi_vs_p", "Photon #phi vs p; Photon Momentum p [GeV]; #Phi [deg]", 
                        50, 0, 4.5, 200, -180, 180);
    TH2D* h_m2_epX_vs_px_miss = new TH2D("h_m2_epX_vs_px_miss", "M_{epX}^{2} vs #DeltaP_{x}; #Delta P_{x} [GeV]; M_{epX}^{2} [GeV^{2}]",
                        100, -0.3, 0.3, 200, -0.5, 0.5);
    TH2D* h_m2_epX_vs_py_miss = new TH2D("h_m2_epX_vs_py_miss", "M_{epX}^{2} vs #Delta P_{y}; #Delta P_{y} [GeV]; M_{epX}^{2} [GeV^{2}]",
                        100, -0.3, 0.3, 200, -0.5, 0.5);
    TH2D* h_m2_epX_vs_pz_miss = new TH2D("h_m2_epX_vs_pz_miss", "M_{epX}^{2} vs #Delta P_{z}; #Delta P_{z} [GeV]; M_{epX}^{2} [GeV^{2}]",
                        100, -0.5, 0.5, 200, -1, 1);
    TH2D* h_m2_epX_vs_E_miss = new TH2D("h_m2_epX_vs_E_miss", "M_{epX}^{2} vs E_{miss}; E_{miss} [GeV]; M_{epX}^{2} [GeV^{2}]",
                        100, -0.5, 0.5, 200, -1, 1);


    // ─── Fill Histograms from Tree ───────────────────────────
    tree->Draw("event.helicity >> h_helicity");
    //tree->Draw("numPhotons >> h_numPhotons");
    tree->Draw("eppi0.m_gg >> h_mgg");
    tree->Draw("dis.Q2 >> h_Q2");
    tree->Draw("dis.W >> h_W");
    tree->Draw("eppi0.t >> h_t");
    tree->Draw("fmod((eppi0.trentoPhi * 180.0 / TMath::Pi()) + 360.0, 360.0) >> h_phiT");
    tree->Draw("eppi0.m2_epX >> h_m2_epX");
    tree->Draw("eppi0.m2_epi0X >> h_m2_epi0X");
    tree->Draw("eppi0.px_miss >> h_px_miss");
    tree->Draw("eppi0.py_miss >> h_py_miss");
    tree->Draw("eppi0.pz_miss >> h_pz_miss");
    tree->Draw("eppi0.E_miss >> h_E_miss");
    tree->Draw("eppi0.pi0_deltaPhi * 180/TMath::Pi() >> h_deltaPhi");

    tree->Draw("dis.Q2:dis.Xb >> h_Q2_vs_Xb", "", "COLZ");
    tree->Draw("eppi0.t: fmod((eppi0.trentoPhi * 180.0 / TMath::Pi()) + 360.0, 360.0) >> h_t_vs_phiT", "", "COLZ");
    tree->Draw("e.theta * 180.0/TMath::Pi() : e.p >> h_e_theta_vs_p", "", "COLZ");
    tree->Draw("e.phi * 180.0/TMath::Pi() : e.p >> h_e_phi_vs_p", "", "COLZ");
    tree->Draw("p.theta * 180.0/TMath::Pi() : p.p >> h_p_theta_vs_p", "", "COLZ");
    tree->Draw("p.phi * 180.0/TMath::Pi() : p.p >> h_p_phi_vs_p", "", "COLZ");
    tree->Draw("g.theta * 180.0/TMath::Pi() : g.p >> h_g_theta_vs_p", "", "COLZ");
    tree->Draw("g.phi * 180.0/TMath::Pi() : g.p >> h_g_phi_vs_p", "", "COLZ");
    tree->Draw("eppi0.m2_epX : eppi0.px_miss >> h_m2_epX_vs_px_miss", "", "COLZ");
    tree->Draw("eppi0.m2_epX : eppi0.py_miss >> h_m2_epX_vs_py_miss", "", "COLZ");
    tree->Draw("eppi0.m2_epX : eppi0.pz_miss >> h_m2_epX_vs_pz_miss", "", "COLZ");
    tree->Draw("eppi0.m2_epX : eppi0.E_miss  >> h_m2_epX_vs_E_miss", "", "COLZ");

    // ─── Fits ──────────────────────────────────
    TF1* mggFit = new TF1("mggFit", "gaus(0) + pol1(3)", 0.098, 0.17);
    Double_t peakPos_mgg = h_mgg->GetBinCenter(h_mgg->GetMaximumBin());
    mggFit->SetParameters(h_mgg->GetMaximum(), peakPos_mgg, 0.01, 10, -1);
    h_mgg->Fit(mggFit, "R");

    TF1* px_miss_fit = new TF1("px_miss_fit", "gaus(0) + pol2(3)", -0.2, 0.2);
    Double_t peakPos_px_miss = h_px_miss->GetBinCenter(h_px_miss->GetMaximumBin());
    px_miss_fit->SetParameters(h_px_miss->GetMaximum(), peakPos_px_miss, 0.05, 1, 1, 1);
    h_px_miss->Fit(px_miss_fit, "R");

    TF1* py_miss_fit = new TF1("py_miss_fit", "gaus(0) + pol2(3)", -0.2, 0.2);
    Double_t peakPos_py_miss = h_py_miss->GetBinCenter(h_py_miss->GetMaximumBin());
    py_miss_fit->SetParameters(h_py_miss->GetMaximum(), peakPos_py_miss, 0.05, 1, 1, 1);
    h_py_miss->Fit(py_miss_fit, "R");

    TF1* pz_miss_fit = new TF1("pz_miss_fit", "gaus(0) + pol2(3)", -0.2, 0.4);
    Double_t peakPos_pz_miss = h_pz_miss->GetBinCenter(h_pz_miss->GetMaximumBin());
    pz_miss_fit->SetParameters(h_pz_miss->GetMaximum(), peakPos_pz_miss, -0.2, 1, 1, 1);
    h_pz_miss->Fit(pz_miss_fit, "R");

    TF1* deltaPhi_fit = new TF1("deltaPhi_fit", "gaus(0) + pol2(3)", -2, 2);
    Double_t peakPos_deltaPhi = h_deltaPhi->GetBinCenter(h_deltaPhi->GetMaximumBin());
    deltaPhi_fit->SetParameters(h_deltaPhi->GetMaximum(), peakPos_deltaPhi, 1, 1, 1, 1);
    h_deltaPhi->Fit(deltaPhi_fit, "R");


    // ─── Create Canvases ─────────────────────────────────────
    TCanvas* c_mgg = new TCanvas("c_mgg", "M_{#gamma#gamma}", 800, 600);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    h_mgg->Draw();
    mggFit->Draw("same");
    c_mgg->Update();

    TCanvas* c_px_miss = new TCanvas("c_px_miss", "#Delta P_{x}", 800, 600);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    h_px_miss->Draw();
    px_miss_fit->Draw("same");
    c_px_miss->Update();

    TCanvas* c_py_miss = new TCanvas("c_py_miss", "#Delta P_{y}", 800, 600);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    h_py_miss->Draw();
    py_miss_fit->Draw("same");
    c_py_miss->Update();

    TCanvas* c_pz_miss = new TCanvas("c_pz_miss", "#Delta P_{z}", 800, 600);
    //gStyle->SetOptStat(1110);
    //gStyle->SetOptFit(1111);
    h_pz_miss->Draw();
    pz_miss_fit->Draw("same");
    c_pz_miss->Update();

    TCanvas* c_deltaPhi = new TCanvas("c_deltaPhi", "#Delta #Phi", 800, 600);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    h_deltaPhi->Draw();
    deltaPhi_fit->Draw("same");
    c_deltaPhi->Update();


    // ─── Save to Output File ─────────────────────────────────
    TFile* outFile = new TFile(outFilePath, "RECREATE");
    h_helicity->Write();
    //h_numPhotons->Write();       
    h_Q2->Write();        
    h_W->Write();
    h_t->Write();
    h_phiT->Write();   
    h_m2_epX->Write();
    h_m2_epi0X->Write();
    h_py_miss->Write();
    h_pz_miss->Write();
    h_E_miss->Write();
    h_deltaPhi->Write();

    h_Q2_vs_Xb->Write();  
    h_t_vs_phiT->Write();
    h_e_theta_vs_p->Write();  
    h_e_phi_vs_p->Write();
    h_p_theta_vs_p->Write();  
    h_p_phi_vs_p->Write();
    h_g_theta_vs_p->Write();
    h_g_phi_vs_p->Write();
    h_m2_epX_vs_px_miss->Write();
    h_m2_epX_vs_py_miss->Write();
    h_m2_epX_vs_pz_miss->Write();
    h_m2_epX_vs_E_miss->Write();

    c_mgg->Write();       
    c_px_miss->Write();
    c_py_miss->Write();
    c_pz_miss->Write();
    c_deltaPhi->Write();

    outFile->Close();

    std::cout << "Histograms and canvases saved to " << outFilePath << std::endl;
}
