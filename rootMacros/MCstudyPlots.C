// To run from command line, type << root -l -b -q 'rootMacros/MCstudyPlots.C("output/test.root")' >>

void MCstudyPlots(const char* inputFilePath = "input.root", const char* outFilePath = "MCstudyPlots.root") {
    // ─── Open Input File and Tree ─────────────────────────────
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

    // ─── Create Histograms ───────────────────────────────────
    // TH1[TYPE]* [h_name] = new TH1[TYPE]("[h_name]", "[TITLE]; [Y TITLE]", [BIN #], [MIN_VAL], [MAX_VAL]);
    // TH2[TYPE]* [h_name] = new TH2[TYPE]("[h_name]", "[y vs x]; [X TITLE] [UNITS]; [Y TITLE] [UNITS]", 
    //                       [X BIN #], [X MIN_VAL], [X MAX_VAL], [Y BIN #], [Y MIN_VAL], [Y MAX_VAL]);


    // INCLUSIVE INFO:
    TH1D* h_pid = new TH1D("h_pid", "Reconstructed particle PID; PID; Counts", 2500, -2500, 2500);

    // DIS INFO:
    TH1D* h_Q2 = new TH1D("h_Q2", "Q^{2}; Q^{2} [GeV^2]; Counts", 100, 0, 6.5);
    TH1D* h_nu = new TH1D("h_nu", "#nu; #nu [GeV]; Counts", 100, 0, 6.5);
    TH1D* h_Xb = new TH1D("h_Xb", "x_{B}; x_{B}; Counts", 100, 0, 1);
    TH1D* h_y  = new TH1D("h_y", "y; y; Counts", 100, 0, 1);
    TH1D* h_W  = new TH1D("h_W", "W; W [GeV]; Counts", 100, 0, 6.5);

    TH1D* h_Q2_gen = new TH1D("h_Q2_gen", "Generated Q^{2}; Q^{2} [GeV^2]; Counts", 100, 0, 6.5);
    TH1D* h_nu_gen = new TH1D("h_nu_gen", "Generated #nu; #nu [GeV]; Counts", 100, 0, 6.5);
    TH1D* h_Xb_gen = new TH1D("h_Xb_gen", "Generated x_{B}; x_{B}; Counts", 100, 0, 1);
    TH1D* h_y_gen  = new TH1D("h_y_gen",  "Generated y; y; Counts", 100, 0, 1);
    TH1D* h_W_gen  = new TH1D("h_W_gen",  "Generated W; W [GeV]; Counts", 100, 0, 6.5);

    TH2D* h_Q2_vs_Xb = new TH2D("h_Q2_vs_Xb", "Reconstructed Q^{2} vs. x_{B}; x_{B} ; Q^{2} [GeV^{2}]",
                                        100, 0, 1, 100, 0, 6.5);

    TH2D* h_Q2_gen_vs_Xb_gen = new TH2D("h_Q2_gen_vs_Xb_gen", "Generated Q^{2} vs. x_{B}; x_{B} ; Q^{2} [GeV^{2}]",
                                        300, 0, 1, 300, 0, 6.5);

    TH2D* h_p_FDdeltaP_vs_p_p = new TH2D("h_p_FDdeltaP_vs_p_p", "FD Proton #Delta p vs. p; p_{p} [GeV]; #Delta p [GeV]",
                                        100, 0, 4, 100, -0.2, 0.2);

    TH2D* h_p_CDdeltaP_vs_p_p = new TH2D("h_p_CDdeltaP_vs_p_p", "CD Proton #Delta p vs. p; p_{p} [GeV]; #Delta p [GeV]",
                                        100, 0, 3, 100, -.1, .1);

    TH2D* h_deltaQ2_vs_Q2 = new TH2D("h_deltaQ2_vs_Q2", "#Delta Q^{2} vs. Q^{2}; Q^{2} [GeV^{2}]; #Delta Q^{2} [GeV^{2}]",
                                        100, 0, 6, 100, -.4, .4);
    
    TH2D* h_deltaXb_vs_Xb = new TH2D("h_deltaXb_vs_Xb", "#Delta x_{B} vs. x_{B}; x_{B}; #Delta x_{B}",
                                        100, 0, 1, 100, -.1, .1);

    TH2D* h_deltat_vs_t = new TH2D("h_deltat_vs_t", "#Delta t vs. t; t [GeV^{2}]; #Delta t [GeV^{2}]",
                                        100, 0, 3, 100, -.4, .4);



    // ─── Fill Histograms from Tree ───────────────────────────
    // tree->Draw("[branch var] >> [h_name]");

    tree->Draw("Q2 >> h_Q2");
    tree->Draw("nu >> h_nu");
    tree->Draw("Xb >> h_Xb");
    tree->Draw("y  >> h_y");
    tree->Draw("W  >> h_W");

    tree->Draw("Q2_gen >> h_Q2_gen");
    tree->Draw("nu_gen >> h_nu_gen");
    tree->Draw("Xb_gen >> h_Xb_gen");
    tree->Draw("y_gen  >> h_y_gen");
    tree->Draw("W_gen  >> h_W_gen");

    tree->Draw("Q2 : Xb >> h_Q2_vs_Xb", "", "COLZ");
    tree->Draw("Q2_gen : Xb_gen >> h_Q2_gen_vs_Xb_gen", "", "COLZ");
    tree->Draw("Q2_gen : Xb_gen >> h_Q2_gen_vs_Xb_gen", "W_gen > 2", "COLZ");

    tree->Draw("p_pgen - p_p : p_p >> h_p_FDdeltaP_vs_p_p", "detPro == 1 && pid == 2212", "COLZ");
    tree->Draw("p_pgen - p_p : p_p >> h_p_CDdeltaP_vs_p_p", "detPro == 2 && pid == 2212", "COLZ");

    tree->Draw("Q2_gen - Q2 : Q2 >> h_deltaQ2_vs_Q2", "", "COLZ");
    tree->Draw("Xb_gen - Xb : Xb >> h_deltaXb_vs_Xb", "", "COLZ");

    tree->Draw("t_gen - t : t >> h_deltat_vs_t", "", "COLZ");
    

    // ─── Fits ──────────────────────────────────
    // TF1* [FIT NAME] = new TF1("[FIT TITLE]", "[FIT FUNC.]", [MIN X VAL], [MAX X VAL]);
    // [FIT NAME]->SetParameters({params});
    // [h_name]->Fit([FIT NAME], "R");


    // ─── Create Canvases ─────────────────────────────────────
    // TCanvas* c = new TCanvas("c", "[TITLE]", 800, 600);
    // gStyle->SetOptStat(1110);
    // gStyle->SetOptFit(1111);
    // h_name->Draw();
    // fit->Draw("same");
    // c->Update();


    // ─── Save to Output File ─────────────────────────────────
    TFile* outFile = new TFile(outFilePath, "RECREATE");

    h_Q2->Write();
    h_nu->Write();
    h_Xb->Write();
    h_y->Write();
    h_W->Write();

    h_Q2_gen->Write();
    h_nu_gen->Write();
    h_Xb_gen->Write();
    h_y_gen->Write();
    h_W_gen->Write();

    h_Q2_vs_Xb->Write();
    h_Q2_gen_vs_Xb_gen->Write();
    
    h_p_FDdeltaP_vs_p_p->Write();
    h_p_CDdeltaP_vs_p_p->Write();
    h_deltaQ2_vs_Q2->Write();
    h_deltaXb_vs_Xb->Write();
    h_deltat_vs_t->Write();

    // c->Write();

    outFile->Close();

    std::cout << "Histograms and canvases saved to " << outFilePath << std::endl;
}
