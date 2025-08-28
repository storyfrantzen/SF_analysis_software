// To run from command line, type << root -l -b -q 'rootMacros/acceptancePlots.C("output/test.root")' >>

void layerBySector(TTree* tree, const char* layerName, const char* yBranch, 
                    const char* xBranch, const char* secBranch) {
    TCanvas* c = new TCanvas(Form("c_%s", layerName), Form("v vs w - %s", layerName), 2000, 1200);
    c->Divide(3, 2);

    for (int sec = 1; sec <= 6; ++sec) {
        c->cd(sec);
        gPad->SetRightMargin(0.15); // Prevents Z-axis numbers from getting cut off
        gPad->SetGrid();
        gPad->SetLogz();
        gStyle->SetOptStat(0); // Disables the stats box for all histograms

        // Define histogram name for reuse
        TString histName = Form("h_%s_vs_%s_sec%d", yBranch, xBranch, sec);

        // Explicit binning: 300x300 bins in v and w from 0 to 450 cm
        TString drawCmd = Form("%s:%s>>%s(450,0,450,450,0,450)", yBranch, xBranch, histName.Data());

        // Cut: sector match AND forward detector electrons only
        TString cut = Form("(%s == %d)", secBranch, sec);

        // Draw and customize
        tree->Draw(drawCmd, cut, "colz");

        char xLabel[4]; // 3 characters + null terminator
        char yLabel[4];
        strncpy(xLabel, xBranch, 3);
        strncpy(yLabel, yBranch, 3);
        xLabel[3] = '\0'; // Ensure null-termination
        yLabel[3] = '\0';

        TH2D* hist = (TH2D*)gPad->FindObject(histName);
        if (hist) {
            hist->SetTitle(Form("%s - Sector %d", layerName, sec));
            hist->GetXaxis()->SetTitle(Form("%s [cm]", xLabel));
            hist->GetYaxis()->SetTitle(Form("%s [cm]", yLabel));
        }
    }

    c->Write();
}

Double_t wrap360(Double_t deg) {
    return deg < 0 ? deg + 360.0 : deg;
}

void acceptancePlots(const char* inputFilePath = "input.root", const char* outFilePath = "acceptancePlots.root") {
    // ─── Open Input File and Tree ─────────────────────────────

    gSystem->Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so");
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

    // INCLUSIVE INFO:
    TH1I* h_pid = new TH1I("h_pid", "Reconstructed particle PID; PID; Counts", 2500, -2500, 2500);
    TH1D* h_helicity = new TH1D("h_helicity", "Helicity Distribution; Helicity; Counts", 10, -5, 5);


    // DIS INFO:
    TH1D* h_Q2 = new TH1D("h_Q2", "Q^{2}; Q^{2} [GeV^2]; Counts", 100, -1, 8);
    TH1D* h_nu = new TH1D("h_nu", "#nu; #nu [GeV]; Counts", 200, 0, 6.5);
    TH1D* h_Xb = new TH1D("h_Xb", "x_{B}; x_{B}; Counts", 200, 0, 1);
    TH1D* h_y  = new TH1D("h_y", "y; y; Counts", 200, 0, 1);
    TH1D* h_W  = new TH1D("h_W", "W; W [GeV]; Counts", 200, 0, 6.5);

    // ELECTRON INFO:
    TH1D* h_e_p = new TH1D("h_e_p", "Electron Momentum; p_{e}; Counts", 300, 0, 6.5);
    TH1D* h_e_theta = new TH1D("h_e_theta", "Electron #theta; #theta_{e}; Counts", 120, 0, 40);
    TH1D* h_e_phi = new TH1D("h_e_phi", "Electron #Phi; #Phi_{e}", 360, 0, 360);
    TH1D* h_e_vz = new TH1D("h_e_vz", "Electron z_{vertex}; z [cm]; Counts", 100, -12, 4);
    TH1D* h_e_chi2pid = new TH1D("h_e_chi2pid", "PID 11 #Chi^{2} Distribution; #Chi^{2}; Counts", 100, 0, 5);
    TH1I* h_detEle = new TH1I("h_detEle", "Electron Detector; Det; Counts", 3, 0, 3);
    TH1I* h_e_sector = new TH1I("h_e_sector", "Electron Sector; Sector; Counts", 6, 1, 7);
    TH1D* h_e_E_PCAL = new TH1D("h_e_E_PCAL", "Electron PCAL Energy; E_{PCAL} [GeV]; Counts", 200, 0, 2);
    TH1D* h_e_E_ECIN = new TH1D("h_e_E_ECIN", "Electron ECIN Energy; E_{ECIN} [GeV]; Counts", 200, 0, 1.5);
    TH1D* h_e_E_ECOUT = new TH1D("h_e_E_ECOUT", "Electron ECOUT Energy; E_{ECOUT} [GeV]; Counts", 200, 0, 1);

    TH2D* h_e_yDC1_vs_e_xDC1 = new TH2D("h_e_yDC1_vs_e_xDC1", "Electron y_{DC1} vs. x_{DC1}; x_{DC1} [cm]; y_{DC1} [cm]",
                                        300, -150, 150, 300, -150, 150);
    TH2D* h_e_yDC2_vs_e_xDC2 = new TH2D("h_e_yDC2_vs_e_xDC2", "Electron y_{DC2} vs. x_{DC2}; x_{DC2} [cm]; y_{DC2} [cm]",
                                        400, -200, 200, 400, -200, 200);
    TH2D* h_e_yDC3_vs_e_xDC3 = new TH2D("h_e_yDC3_vs_e_xDC3", "Electron y_{DC3} vs. x_{DC3}; x_{DC3} [cm]; y_{DC3} [cm]",
                                        400, -300, 300, 400, -300, 300);

    TH2D* h_Q2_vs_Xb = new TH2D("h_Q2_vs_Xb", "Q^{2} vs. x_{B}; x_{B}; Q^{2} [GeV^{2}]",
                                        500, 0, 1, 500, 0, 7);

    // PROTON INFO:
    TH1D* h_p_p = new TH1D("h_p_p", "Proton Momentum; p_{p}; Counts", 100, 0, 6.5);
    TH1D* h_p_theta = new TH1D("h_p_theta", "Proton #theta; #theta_{p}; Counts", 200, 0, 100);
    TH1D* h_p_phi = new TH1D("h_p_phi", "Proton #Phi; #Phi_{p}; Counts", 360, 0, 360);
    TH1D* h_p_vz = new TH1D("h_p_vz", "Proton z_{vertex}; z [cm]; Counts", 100, -12, 4);
    TH1D* h_p_chi2pid = new TH1D("h_p_chi2pid", "PID 2212 #Chi^{2} Distribution; #Chi^{2}; Counts", 100, 0, 5);
    TH1I* h_detPro = new TH1I("h_detPro", "Proton Detector; Det; Counts", 3, 0, 3);
    TH1I* h_p_sector = new TH1I("h_p_sector", "Proton Sector; Sector; Counts", 6, 1, 7);
    TH1D* h_p_edge_cvt1 = new TH1D("h_p_edge_cvt1", "Proton CVT1 Edge; Edge [cm]; Counts", 8, -4, 4);
    TH1D* h_p_edge_cvt3 = new TH1D("h_p_edge_cvt3", "Proton CVT3 Edge; Edge [cm]; Counts", 8, -4, 4);
    TH1D* h_p_edge_cvt5 = new TH1D("h_p_edge_cvt5", "Proton CVT5 Edge; Edge [cm]; Counts", 8, -4, 4);
    TH1D* h_p_edge_cvt7 = new TH1D("h_p_edge_cvt7", "Proton CVT7 Edge; Edge [cm]; Counts", 14, -4, 10);
    TH1D* h_p_edge_cvt12 = new TH1D("h_p_edge_cvt12", "Proton CVT12 Edge; Edge [cm]; Counts", 14, -4, 10);
    TH1D* h_p_theta_cvt = new TH1D("h_p_theta_cvt", "Proton #Theta_{CVT}; #Theta_{CVT} [deg]; Counts", 180, 0, 180);
    TH1D* h_p_phi_cvt = new TH1D("h_p_phi_cvt", "Proton #Phi_{CVT}; #Phi_{CVT} [deg]; Counts", 360, 0, 360);

    TH2D* h_p_yDC1_vs_p_xDC1 = new TH2D("h_p_yDC1_vs_p_xDC1", "Proton y_{DC1} vs. x_{DC1}; x_{DC1} [cm]; y_{DC1} [cm]",
                                        400, -200, 200, 400, -200, 200);
    TH2D* h_p_yDC2_vs_p_xDC2 = new TH2D("h_p_yDC2_vs_p_xDC2", "Proton y_{DC2} vs. x_{DC2}; x_{DC2} [cm]; y_{DC2} [cm]",
                                        500, -250, 250, 500, -250, 250);
    TH2D* h_p_yDC3_vs_p_xDC3 = new TH2D("h_p_yDC3_vs_p_xDC3", "Proton y_{DC3} vs. x_{DC3}; x_{DC3} [cm]; y_{DC3} [cm]",
                                        700, -350, 350, 700, -350, 350);

    TH2D* h_p_theta_vs_p_phi = new TH2D("h_p_theta_vs_p_phi", "Proton #theta vs. #Phi; #Phi_{p} [deg]; #theta_{p} [deg]",
                                        360, 0, 360, 140, 0, 140);
    TH2D* h_p_theta_cvt_vs_p_phi_cvt = new TH2D("h_p_theta_cvt_vs_p_phi_cvt", "Proton #theta_{CVT} vs. #Phi_{CVT}; #Phi_{CVT} [deg]; #theta_{CVT} [deg]", 
                                        360, 0, 360, 140, 0, 140);
            

    // PHOTON INFO:
    TH1D* h_g_p = new TH1D("h_g_p", "Photon Momentum; g_{p}; Counts", 100, 0, 6.5);
    TH1D* h_g_theta = new TH1D("h_g_theta", "Photon #theta; #theta_{#gamma}; Counts", 200, 0, 100);
    TH1D* h_g_phi = new TH1D("h_g_phi", "Photon #Phi; #Phi_{p}; Counts", 360, 0, 360);
    TH1D* h_g_chi2pid = new TH1D("h_g_chi2pid", "PID 22 #Chi^{2} Distribution; #Chi^{2}; Counts", 100, 0, 5);
    TH1I* h_detPho = new TH1I("h_detPho", "Photon Detector; Det; Counts", 2, 0, 2);
    TH1I* h_g_sector = new TH1I("h_g_sector", "Photon Sector; Sector; Counts", 6, 1, 7);
    TH1D* h_g_E_PCAL = new TH1D("h_g_E_PCAL", "Photon PCAL Energy; E_{PCAL} [GeV]; Counts", 100, 0, 2);
    TH1D* h_g_E_ECIN = new TH1D("h_g_E_ECIN", "Photon ECIN Energy; E_{ECIN} [GeV]; Counts", 100, 0, 1);
    TH1D* h_g_E_ECOUT = new TH1D("h_g_E_ECOUT", "Photon ECOUT Energy; E_{ECOUT} [GeV]; Counts", 100, 0, 1);

    

    // ─── Fill Histograms from Tree ───────────────────────────
    tree->Draw("event.helicity >> h_helicity");

    tree->Draw("dis.Q2 >> h_Q2");
    tree->Draw("dis.nu >> h_nu");
    tree->Draw("dis.Xb >> h_Xb");
    tree->Draw("dis.y  >> h_y");
    tree->Draw("dis.W  >> h_W");

    tree->Draw("e.p >> h_e_p");
    tree->Draw("e.theta * 180.0/TMath::Pi() >> h_e_theta");
    tree->Draw("wrap360(e.phi * 180.0/TMath::Pi()) >> h_e_phi");
    tree->Draw("e.vz >> h_e_vz");
    tree->Draw("e.chi2pid >> h_e_chi2pid");
    tree->Draw("e.det >> h_detEle", "e.det >= 0");
    tree->Draw("e.sector >> h_e_sector", "e.sector > -1");
    tree->Draw("e.E_PCAL >> h_e_E_PCAL");
    tree->Draw("e.E_ECIN >> h_e_E_ECIN");
    tree->Draw("e.E_ECOUT >> h_e_E_ECOUT");

    tree->Draw("e.yDC1:e.xDC1 >> h_e_yDC1_vs_e_xDC1", "", "COLZ");
    tree->Draw("e.yDC2:e.xDC2 >> h_e_yDC2_vs_e_xDC2", "", "COLZ");
    tree->Draw("e.yDC3:e.xDC3 >> h_e_yDC3_vs_e_xDC3", "", "COLZ");

    tree->Draw("dis.Q2 : dis.Xb >> h_Q2_vs_Xb", "", "COLZ");

    tree->Draw("p.p >> h_p_p");
    tree->Draw("p.theta * 180.0/TMath::Pi() >> h_p_theta");
    tree->Draw("wrap360(p.phi * 180.0/TMath::Pi()) >> h_p_phi");
    tree->Draw("p.vz >> h_p_vz");
    tree->Draw("p.chi2pid >> h_p_chi2pid");
    tree->Draw("p.det >> h_detPro", "p.det >= 0");
    tree->Draw("p.sector >> h_p_sector", "p.sector > -1");
    tree->Draw("p.edge_cvt1 >> h_p_edge_cvt1");
    tree->Draw("p.edge_cvt3 >> h_p_edge_cvt3");
    tree->Draw("p.edge_cvt5 >> h_p_edge_cvt5");
    tree->Draw("p.edge_cvt7 >> h_p_edge_cvt7");
    tree->Draw("p.edge_cvt12 >> h_p_edge_cvt12");
    tree->Draw("p.theta_cvt * 180/TMath::Pi() >> h_p_theta_cvt");
    tree->Draw("wrap360(p.phi_cvt * 180/TMath::Pi()) >> h_p_phi_cvt", "TMath::Finite(p.phi_cvt)");

    tree->Draw("p.yDC1:p.xDC1 >> h_p_yDC1_vs_p_xDC1", "", "COLZ");
    tree->Draw("p.yDC2:p.xDC2 >> h_p_yDC2_vs_p_xDC2", "", "COLZ");
    tree->Draw("p.yDC3:p.xDC3 >> h_p_yDC3_vs_p_xDC3", "", "COLZ");
    tree->Draw("p.theta * 180.0/TMath::Pi() : wrap360(p.phi * 180.0/TMath::Pi()) >> h_p_theta_vs_p_phi");
    tree->Draw("p.theta_cvt * 180.0/TMath::Pi() : wrap360(p.phi_cvt * 180.0/TMath::Pi()) >> h_p_theta_cvt_vs_p_phi_cvt");

    tree->Draw("g.p >> h_g_p");
    tree->Draw("g.theta * 180.0/TMath::Pi() >> h_g_theta");
    tree->Draw("wrap360(g.phi * 180.0/TMath::Pi()) >> h_g_phi");
    tree->Draw("g.chi2pid >> h_g_chi2pid");
    tree->Draw("g.det >> h_detPho", "g.det >= 0");
    tree->Draw("g.sector >> h_g_sector", "g.sector > -1");
    tree->Draw("g.E_PCAL >> h_g_E_PCAL");
    tree->Draw("g.E_ECIN >> h_g_E_ECIN");
    tree->Draw("g.E_ECOUT >> h_g_E_ECOUT");

    
    // ─── Save to Output File ─────────────────────────────────
    TFile* outFile = new TFile(outFilePath, "RECREATE");

    //[h_name]->Write();
    h_pid->Write();
    h_helicity->Write();

    h_Q2->Write();
    h_nu->Write();
    h_Xb->Write();
    h_y->Write();
    h_W->Write();

    h_e_p->Write();
    h_e_theta->Write();
    h_e_phi->Write();
    h_e_vz->Write();
    h_e_chi2pid->Write();
    h_detEle->Write();
    h_e_sector->Write();
    h_e_E_PCAL->Write();
    h_e_E_ECIN->Write();
    h_e_E_ECOUT->Write();

    h_e_yDC1_vs_e_xDC1->Write();
    h_e_yDC2_vs_e_xDC2->Write();
    h_e_yDC3_vs_e_xDC3->Write();

    h_Q2_vs_Xb->Write();

    layerBySector(tree, "PCAL", "e.vPCAL", "e.wPCAL", "e.sector");
    layerBySector(tree, "ECIN", "e.vECIN", "e.wECIN", "e.sector");
    layerBySector(tree, "ECOUT", "e.vECOUT", "e.wECOUT", "e.sector");
    layerBySector(tree, "vu_ECOUT", "e.vECOUT", "e.uECOUT", "e.sector");

    h_p_p->Write();
    h_p_theta->Write();
    h_p_phi->Write();
    h_p_vz->Write();
    h_p_chi2pid->Write();
    h_detPro->Write();
    h_p_sector->Write();
    h_p_edge_cvt1->Write();
    h_p_edge_cvt3->Write();
    h_p_edge_cvt5->Write();
    h_p_edge_cvt7->Write();
    h_p_edge_cvt12->Write();
    h_p_theta_cvt->Write();
    h_p_phi_cvt->Write();

    h_p_yDC1_vs_p_xDC1->Write();
    h_p_yDC2_vs_p_xDC2->Write();
    h_p_yDC3_vs_p_xDC3->Write();
    h_p_theta_vs_p_phi->Write();
    h_p_theta_cvt_vs_p_phi_cvt->Write();

    h_g_p->Write();
    h_g_theta->Write();
    h_g_phi->Write();
    h_g_chi2pid->Write();
    h_detPho->Write();
    h_g_sector->Write();
    h_g_E_PCAL->Write();
    h_g_E_ECIN->Write();
    h_g_E_ECOUT->Write();

    layerBySector(tree, "gPCAL", "g.vPCAL", "g.wPCAL", "g.sector");
    layerBySector(tree, "gECIN", "g.vECIN", "g.wECIN", "g.sector");
    layerBySector(tree, "gECOUT", "g.vECOUT", "g.wECOUT", "g.sector");
    layerBySector(tree, "vu_gECOUT", "g.vECOUT", "g.uECOUT", "g.sector");

    // c->Write();

    outFile->Close();

    std::cout << "Histograms and canvases saved to " << outFilePath << std::endl;
}
