    // To run from command line, type << root -l -b -q 'rootMacros/analyzeBSA.C("output/test.root")' >>

    #include <cmath>

    void analyzeBSA(const char* inputFilePath = "/work/clas12/storyf/SF_analysis_software/output/test.root") {
        TFile *f = TFile::Open(inputFilePath);
        if (!f || f->IsZombie()) {
            std::cerr << "Could not open input file: " << inputFilePath << std::endl;
            return;
        }
        TTree *tree = (TTree*)f->Get("Events");
        if (!tree) {
            std::cerr << "TTree 'Events' not found in file: " << inputFilePath << std::endl;
            return;
        }

        double Q2, Xb, t, phiT;
        int helicity;
        // Link tree branches to post-process variables, e.g., tree->SetBranchAddress("branchName", &branchVar);
        tree->SetBranchAddress("Q2", &Q2);
        tree->SetBranchAddress("Xb", &Xb);
        tree->SetBranchAddress("t",  &t);
        tree->SetBranchAddress("trentoPhi", &phiT);
        tree->SetBranchAddress("helicity", &helicity);

        // --- Step 1: Adaptive binning for Q2, Xb, t ---
        auto getAdaptiveEdges = [&](const char* var, int nBins, double min, double max) {
            TH1D* h = new TH1D(Form("h_%s", var), "", 500, min, max);
            tree->Draw(Form("%s >> h_%s", var, var), "", "goff");
            std::vector<double> edges;
            int total = h->GetEntries();
            int step = total / nBins;
            //std::cout << "step = " << step << std::endl;
            int cumulative = 0;
            edges.push_back(h->GetXaxis()->GetXmin());

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                cumulative += h->GetBinContent(i);
                if (cumulative >= step * edges.size() && edges.back() < h->GetBinLowEdge(i+1)) {
                    edges.push_back(h->GetBinLowEdge(i+1));
                }
            }

            // Ensure upper edge included
            if (edges.back() < h->GetXaxis()->GetXmax())
                edges.push_back(h->GetXaxis()->GetXmax());

            delete h;
            return edges;
        };

        // get vector lists of optimal edges for each kinematic quantity, but uniform phi binning
        std::vector<double> Q2_edges     = getAdaptiveEdges("Q2", 4, 1, 4);
        std::vector<double> Xb_edges     = getAdaptiveEdges("Xb", 4, 0.1, 0.6);
        std::vector<double> t_edges      = getAdaptiveEdges("t",  1, 0, 3);
        std::vector<double> phi_edges;
        for (int i = 0; i <= 12; ++i) phi_edges.push_back(i * 30.0); // Uniform phi bins

        const int numQ2  = Q2_edges.size() - 1;
        const int numXb  = Xb_edges.size() - 1;
        const int numt   = t_edges.size()  - 1;
        const int numphi = phi_edges.size() - 1;
        double beamPol = 0.8;

        std::cout << "numQ2 = " << numQ2 << ", numXb = " << numXb << ", numt = " << numt << ", numphi = " << numphi << std::endl;

        auto getBinIndex = [](double value, const std::vector<double>& edges) -> int {
            for (size_t i = 0; i < edges.size() - 1; ++i) {
                if (value >= edges[i] && value < edges[i + 1]) return i;
            }
            return -1;
        };

        // 1D arrays for phi-differential BSA
        std::vector<int> Np_int(numphi, 0);
        std::vector<int> Nm_int(numphi, 0);

        // 4D arrays for fully differential counts
        std::vector<std::vector<std::vector<std::vector<int>>>> Np(numQ2,
            std::vector<std::vector<std::vector<int>>>(numXb,
                std::vector<std::vector<int>>(numt,
                    std::vector<int>(numphi, 0))));

        auto Nm = Np;

        // 4D arrays for <Q2>, <Xb>, <t> in given kinematic bin
        std::vector<std::vector<std::vector<std::vector<double>>>> sum_Q2(numQ2, std::vector<std::vector<std::vector<double>>>(numXb, std::vector<std::vector<double>>(numt, std::vector<double>(numphi, 0.0))));
        std::vector<std::vector<std::vector<std::vector<double>>>> sum_Xb = sum_Q2;
        std::vector<std::vector<std::vector<std::vector<double>>>> sum_t  = sum_Q2;

        // Loop over events
        Long64_t numEntries = tree->GetEntries();
        for (Long64_t i = 0; i < numEntries; ++i) {
            tree->GetEntry(i);

            double phi_deg = fmod(phiT * 180.0 / TMath::Pi() + 360.0, 360.0);

            int iQ2  = getBinIndex(Q2, Q2_edges);
            int iXb  = getBinIndex(Xb, Xb_edges);
            int it   = getBinIndex(t, t_edges); 
            int iphi = getBinIndex(phi_deg, phi_edges);

            if (iQ2 < 0 || iXb < 0 || it < 0 || iphi < 0) continue;

            if (helicity == +1 || helicity == -1) {
                sum_Q2[iQ2][iXb][it][iphi] += Q2;
                sum_Xb[iQ2][iXb][it][iphi] += Xb;
                sum_t[iQ2][iXb][it][iphi]  += t;
            }

            if (helicity == +1) {
                Np_int[iphi]++;
                Np[iQ2][iXb][it][iphi]++;
            }
            else if (helicity == -1) {
                Nm_int[iphi]++;
                Nm[iQ2][iXb][it][iphi]++;
            }
        }

        for (int i = 0; i < numQ2; ++i)
            for (int j = 0; j < numXb; ++j) 
                for (int k = 0; k < numt; ++k)
                    for (int l = 0; l < numphi; ++l) {
                        int total = Np[i][j][k][l] + Nm[i][j][k][l];
                        if (total > 0)
                        std::cout << Form("Bin Q2[%d] Xb[%d] t[%d] phi[%d]: N+ = %d, N- = %d\n", i, j, k, l, Np[i][j][k][l], Nm[i][j][k][l]);
                    }


        // Output results
        TH1D* h_phi_asym[numQ2][numXb][numt];
        for (int i = 0; i < numQ2; ++i) 
            for (int j = 0; j < numXb; ++j) 
                for (int k = 0; k < numt; ++k) {
                    TString name = Form("Q2%d_Xb%d_t%d", i, j, k);
                    h_phi_asym[i][j][k] = new TH1D(name, name, numphi, &phi_edges[0]);
                    
                    for (int l = 0; l < numphi; ++l) {
                        int Nplus  = Np[i][j][k][l];
                        int Nminus = Nm[i][j][k][l];
                        int Nsum   = Nplus + Nminus;

                        if (Nsum > 0) {
                            double ALU = (Nplus - Nminus) / (double)Nsum / beamPol;
                            h_phi_asym[i][j][k]->SetBinContent(l + 1, ALU);
                        }
                    }
                }
            
        

        TH1D* h_phi_asym_int = new TH1D("A_LU_int", "A_{LU} vs #phi (all kinematics); #phi [deg]; A_{LU}", numphi, &phi_edges[0]);
        TF1* fit_asym_int = new TF1("fit_asym_int", "[0] * sin([1] * x * TMath::DegToRad() + [2])", 0, 360);
        fit_asym_int->SetParameters(0.2, 1.0, 0);

        for (int l = 0; l < numphi; ++l) {
            int Nplus  = Np_int[l];
            int Nminus = Nm_int[l];
            int Nsum   = Nplus + Nminus;

            // Not sure if min Nsum needed here
            if (Nsum > 0) {
                double ALU = (Nplus - Nminus) / (double)Nsum / beamPol;

                double Nplus_d  = static_cast<double>(Nplus);
                double Nminus_d = static_cast<double>(Nminus);
                double Nsum_d   = Nplus_d + Nminus_d;

                double error = (2.0 / beamPol) * sqrt((Nplus_d * Nminus_d) / pow(Nsum_d, 3));

                std::cout << "Bin " << l << " has ALU = " << ALU << " and error = " << error << std::endl;

                h_phi_asym_int->SetBinContent(l + 1, ALU);
                h_phi_asym_int->SetBinError(l + 1, error);

            }
        }

        h_phi_asym_int->Fit(fit_asym_int, "R");

        // Optional: save to file
        TFile fout("bsa_output.root", "RECREATE");

        // write separate BSA for each kinematic bin, d4sigma/dQ2/dXb/dt/dphi : 
        for (int i = 0; i < numQ2; ++i)
            for (int j = 0; j < numXb; ++j)
                for (int k = 0; k < numt; ++k)
                    h_phi_asym[i][j][k]->Write();

        // write the integrated BSA, dsigma/dphi :
        h_phi_asym_int->Draw("PE1");
        h_phi_asym_int->Write();

        // construct and populate canvas for 2D binning of d2sigma/dQ2/dXb :
        TCanvas* c_2D = new TCanvas("c_2D", "2D BSA", 3000, 2000);
        c_2D->cd();

        double dx = 1.0 / (numXb + 0.1);
        double dy = 1.0 / (numQ2 + 0.7);

        for (int i = 0; i < numQ2; ++i) {
            int iq_flipped = numQ2 - 1 - i;
            for (int j = 0; j < numXb; ++j) {
                bool hasEntries = false;
                for (int k = 0; k < numt; ++k) {
                    if (h_phi_asym[i][j][k] && h_phi_asym[i][j][k]->GetEntries() > 5) {
                        hasEntries = true;
                        break;
                    }
                }
                if (!hasEntries) continue;

                double x0 = 0.05 + j *  dx;
                double y0 = 0.1 + (numQ2 - 1 - iq_flipped) * dy;
                double x1 = x0 + dx ;
                double y1 = y0 + dy;

                TString padName = Form("pad_Q2_%d_xB_%d", i, j);
                TPad* pad = new TPad(padName, padName, x0, y0, x1, y1);
                pad->SetMargin(0.2, 0.05, 0.2, 0.05);
                pad->Draw();
                pad->cd();

                // Draw the first t-bin only (can be adapted)
                // Customize plot appearance for clarity
                TH1D* h = h_phi_asym[i][j][0];
                h->SetStats(0);
                h->SetTitle("");
                h->GetXaxis()->SetTitle("#phi [deg]");
                h->GetYaxis()->SetTitle("A_{LU}");
                h->GetXaxis()->CenterTitle();
                h->GetYaxis()->CenterTitle();
                h->GetXaxis()->SetTitleSize(0.08);
                h->GetYaxis()->SetTitleSize(0.08);
                h->GetXaxis()->SetLabelSize(0.07);
                h->GetYaxis()->SetLabelSize(0.07);
                h->SetMarkerStyle(20);
                h->SetMarkerSize(.5); // Increase marker size for visibility
                h->Draw("P E1");

                // Label each plot with Q2/xB info
                TLatex label;
                label.SetNDC(true);
                label.SetTextAlign(13);
                label.SetTextSize(0.07);
                label.DrawLatex(0.1, 0.85, Form("Q^{2}_{bin} %d", i));
                label.DrawLatex(0.1, 0.75, Form("x_{B}_{bin} %d", j));

                c_2D->cd();
            }
        }
        c_2D->cd();
        c_2D->Update();
        c_2D->Write();          

        fout.Close();
    }
