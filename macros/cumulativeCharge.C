// cumulativeCharge.C
void cumulativeCharge() {
    TFile *f = TFile::Open("output/RGK_FID_SKIM18_eppi0REC_6.535_tor1_1028_1601.root");
    TTree *t = (TTree*)f->Get("Events");

    Long64_t n = t->GetEntries();
    std::vector<double> cumCharge(n);
    std::vector<double> cumEvents(n);

    double chargeSum = 0;
    for (Long64_t i=0; i<n; ++i) {
        double q;
        t->GetEntry(i);
        t->SetBranchAddress("event.charge", &q); 
        chargeSum += q;
        cumCharge[i] = chargeSum;
        cumEvents[i] = i+1;
    }

    TGraph *g = new TGraph(n, &cumCharge[0], &cumEvents[0]);
    g->SetTitle("Cumulative Events vs. Cumulative Charge;Integrated Charge;Number of Events");
    g->SetLineColor(kBlue);
    g->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1","Cumulative Plot",800,600);
    g->Draw("AP"); // A=axis, L=line
    c1->SaveAs("cumulativeCharge.png");
}
