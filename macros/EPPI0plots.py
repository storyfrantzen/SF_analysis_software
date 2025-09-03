#type: ignore
import ROOT
import argparse
import sys

def check_branches(tree, branches):
    missing = [b for b in branches if not tree.GetBranch(b)]
    if missing:
        print("\nERROR: Missing required branches in the TTree 'Events':")
        for mb in missing:
            print(f"  - {mb}")
        print("\nPlease ensure the input ROOT file contains all required branches.")
        print("Required branches:")
        for b in branches:
            print(f"  {b}")
        return False
    return True


def EPPI0plots(inputFilePath="input.root", outFilePath="EPPI0plots.root"):
    ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

    # Required branches
    required_branches = [
        "event.helicity",
        "dis.Q2","dis.nu","dis.Xb","dis.y","dis.W",
        "e.p","e.theta","e.phi","e.vz","e.chi2pid","e.det","e.sector",
        "e.E_PCAL","e.E_ECIN","e.E_ECOUT",
        "e.xDC1","e.yDC1","e.xDC2","e.yDC2","e.xDC3","e.yDC3",
        "p.p","p.theta","p.phi","p.vz","p.chi2pid","p.det","p.sector",
        "p.edge_cvt1","p.edge_cvt3","p.edge_cvt5","p.edge_cvt7","p.edge_cvt12",
        "p.theta_cvt","p.phi_cvt",
        "p.xDC1","p.yDC1","p.xDC2","p.yDC2","p.xDC3","p.yDC3",
        "g.p","g.theta","g.phi","g.chi2pid","g.det","g.sector",
        "g.E_PCAL","g.E_ECIN","g.E_ECOUT",
        "e.vPCAL","e.wPCAL","e.vECIN","e.wECIN","e.vECOUT","e.wECOUT","e.uECOUT",
        "g.vPCAL","g.wPCAL","g.vECIN","g.wECIN","g.vECOUT","g.wECOUT","g.uECOUT"
    ]

    f = ROOT.TFile.Open(inputFilePath)
    if not f or f.IsZombie():
        print(f"Could not open file: {inputFilePath}")
        return

    tree = f.Get("Events")
    if not tree:
        print("Could not find tree 'Events' in file.")
        return

    check_branches(tree, required_branches)

    # ─── Histograms ─────────────────────────────
    histos1D = {
        "h_helicity": ROOT.TH1I("h_helicity", "Helicity; Counts", 10, -5, 5),
        "h_mgg_pFD": ROOT.TH1D("h_mgg_pFD", "pFD M_{#gamma#gamma}; M_{#gamma#gamma} [GeV]; Counts", 200, 0.098, 0.17),
        "h_mgg_pCD": ROOT.TH1D("h_mgg_pCD", "pCD M_{#gamma#gamma}; M_{#gamma#gamma} [GeV]; Counts", 200, 0.098, 0.17),
        "h_Q2": ROOT.TH1D("h_Q2", "Q^{2} Distribution; Q^{2} [GeV^{2}]; Counts", 250, 0, 10),
        "h_W": ROOT.TH1D("h_W", "W Distribution; W [GeV]; Counts", 100, 2, 3.5),
        "h_t": ROOT.TH1D("h_t", "-t Distribution; -t [GeV^{2}]; Counts", 100, 0, 4),
        "h_phiT": ROOT.TH1D("h_phiT", "Trento #Phi Distribution; #Phi_{T} [deg]; Counts", 180, 0, 360),
        "h_theta_e_g1": ROOT.TH1D("h_theta_e_g1", "#theta_{e#gamma_{1}}; #theta [deg]; Counts", 90, 0, 90),
        "h_theta_e_g2": ROOT.TH1D("h_theta_e_g2", "#theta_{e#gamma_{2}}; #theta [deg]; Counts", 90, 0, 90),
        "h_theta_g1_g2": ROOT.TH1D("h_theta_g1_g2", "#theta_{#gamma_{1}#gamma_{2}}; #theta [deg]; Counts", 90, 0, 90),
        "h_m2_miss_pFD": ROOT.TH1D("h_m2_miss_pFD", "pFD M_{X}^{2}; M_{X}^{2} [GeV^{2}]; Counts", 200, -1, 4),
        "h_m2_miss_pCD": ROOT.TH1D("h_m2_miss_pCD", "pCD M_{X}^{2}; M_{X}^{2} [GeV^{2}]; Counts", 200, -1, 4),
        "h_m2_epX_pFD": ROOT.TH1D("h_m2_epX_pFD", "pFD M_{epX}^{2}; M_{epX}^{2} [GeV^{2}]; Counts", 200, -1, 4),
        "h_m2_epX_pCD": ROOT.TH1D("h_m2_epX_pCD", "pCD M_{epX}^{2}; M_{epX}^{2} [GeV^{2}]; Counts", 200, -1, 4),
        "h_m2_epi0X_pFD": ROOT.TH1D("h_m2_epi0X_pFD", "pFD M_{e#piX}^{2}; M_{e#piX}^{2} [GeV^{2}]; Counts", 200, -1, 4),
        "h_m2_epi0X_pCD": ROOT.TH1D("h_m2_epi0X_pCD", "pCD M_{e#piX}^{2}; M_{e#piX}^{2} [GeV^{2}]; Counts", 200, -1, 4),
        "h_px_miss": ROOT.TH1D("h_px_miss", "Missing Momentum #Delta P_{x}; #Delta P_{x} [GeV]; Counts", 200, -0.8, 0.8),
        "h_py_miss": ROOT.TH1D("h_py_miss", "Missing Momentum #Delta P_{y}; #Delta P_{y} [GeV]; Counts", 200, -0.8, 0.8),
        "h_pz_miss": ROOT.TH1D("h_pz_miss", "Missing Momentum #Delta P_{z}; #Delta P_{z} [GeV]; Counts", 200, -0.3, 0.6),
        "h_E_miss_pFD": ROOT.TH1D("h_E_miss_pFD", "pFD E_{miss}; E_{miss} [GeV]; Counts", 200, -0.8, 1.2),
        "h_E_miss_pCD": ROOT.TH1D("h_E_miss_pCD", "pCD E_{miss}; E_{miss} [GeV]; Counts", 200, -0.8, 1.2),
        "h_deltaPhi": ROOT.TH1D("h_deltaPhi", "#Delta #Phi; #Delta #Phi [deg]; Counts", 200, -10, 10)
    }

    # 2D histograms
    histos2D = {
        "h_Q2_vs_Xb": ROOT.TH2D("h_Q2_vs_Xb", "Q^{2} vs x_{B}; x_{B}; Q^{2} [GeV^{2}]", 50, 0, 0.7, 50, 0.5, 7),
        "h_Q2_vs_Xb_pFD": ROOT.TH2D("h_Q2_vs_Xb_pFD", "pFD Q^{2} vs x_{B}; x_{B}; Q^{2} [GeV^{2}]", 50, 0, 0.7, 50, 0.5, 7),
        "h_Q2_vs_Xb_pCD": ROOT.TH2D("h_Q2_vs_Xb_pCD", "pCD Q^{2} vs x_{B}; x_{B}; Q^{2} [GeV^{2}]", 50, 0, 0.7, 50, 0.5, 7),
        "h_t_vs_phiT": ROOT.TH2D("h_t_vs_phiT", "-t vs. #Phi_{T}; #Phi_{T} [deg]; -t [GeV^{2}]", 180, 0, 360, 50, 0, 2),
        "h_e_theta_vs_p": ROOT.TH2D("h_e_theta_vs_p", "Electron #theta vs p; Electron Momentum p [GeV]; #theta [deg]", 50, 1, 4.5, 80, 0, 40),
        "h_e_phi_vs_p": ROOT.TH2D("h_e_phi_vs_p", "Electron #phi vs p; Electron Momentum p [GeV]; #Phi [deg]", 50, 1, 4.5, 180, -180, 180),
        "h_p_theta_vs_p": ROOT.TH2D("h_p_theta_vs_p", "Proton #theta vs p; Proton Momentum p [GeV]; #theta [deg]", 50, 0, 4.5, 100, 0, 90),
        "h_p_phi_vs_p": ROOT.TH2D("h_p_phi_vs_p", "Proton #phi vs p; Proton Momentum p [GeV]; #Phi [deg]", 50, 0, 4.5, 100, -180, 180),
        "h_g_theta_vs_p": ROOT.TH2D("h_g_theta_vs_p", "Photon #theta vs p; Photon Momentum p [GeV]; #theta [deg]", 50, 0, 5, 100, 0, 40),
        "h_g_phi_vs_p": ROOT.TH2D("h_g_phi_vs_p", "Photon #phi vs p; Photon Momentum p [GeV]; #Phi [deg]", 50, 0, 5, 180, -180, 180),
        "h_m2_miss_vs_p_theta": ROOT.TH2D("h_m2_miss_vs_p_theta", "M_{X}^{2} vs #theta_{p}; #theta_{p} [deg]; M_{X}^{2} [GeV^{2}]", 100, 0, 70, 100, -1, 2),
        "h_m2_epX_vs_p_theta": ROOT.TH2D("h_m2_epX_vs_p_theta", "M_{epX}^{2} vs #theta_{p}; #theta_{p} [deg]; M_{epX}^{2} [GeV^{2}]", 100, 0, 70, 100, -1, 2),
        "h_m2_epi0X_vs_p_theta": ROOT.TH2D("h_m2_epi0X_vs_p_theta", "M_{e#piX}^{2} vs #theta_{p}; #theta_{p} [deg]; M_{e#piX}^{2} [GeV^{2}]", 100, 0, 70, 100, -1, 3),
        "h_m2_miss_vs_px_miss": ROOT.TH2D("h_m2_miss_vs_px_miss", "M_{X}^{2} vs #DeltaP_{x}; #Delta P_{x} [GeV]; M_{X}^{2} [GeV^{2}]", 100, -0.3, 0.3, 200, -0.5, 0.5),
        "h_m2_miss_vs_py_miss": ROOT.TH2D("h_m2_miss_vs_py_miss", "M_{X}^{2} vs #DeltaP_{y}; #Delta P_{y} [GeV]; M_{X}^{2} [GeV^{2}]", 100, -0.3, 0.3, 200, -0.5, 0.5),
        "h_m2_miss_vs_pz_miss": ROOT.TH2D("h_m2_miss_vs_pz_miss", "M_{X}^{2} vs #Delta P_{z}; #Delta P_{z} [GeV]; M_{X}^{2} [GeV^{2}]", 100, -0.5, 0.5, 200, -1, 1),
        "h_m2_miss_vs_E_miss": ROOT.TH2D("h_m2_miss_vs_E_miss", "M_{X}^{2} vs E_{miss}; E_{miss} [GeV]; M_{X}^{2} [GeV^{2}]", 100, -0.5, 0.5, 200, -1, 1)
    }

    # ─── Fill histograms ─────────────────────────
    # Example: adjust selection cuts as needed (from original macro)
    tree.Draw("event.helicity >> h_helicity")
    tree.Draw("eppi0.m_gg >> h_mgg_pFD", "p.det==1")
    tree.Draw("eppi0.m_gg >> h_mgg_pCD", "p.det==2")
    tree.Draw("dis.Q2 >> h_Q2")
    tree.Draw("dis.W >> h_W")
    tree.Draw("t >> h_t")
    tree.Draw("fmod((eppi0.trentoPhi * 180.0 / TMath::Pi()) + 360.0, 360.0) >> h_phiT")
    tree.Draw("eppi0.theta_e_g1 * 180.0/TMath::Pi() >> h_theta_e_g1")
    tree.Draw("eppi0.theta_e_g2 * 180.0/TMath::Pi() >> h_theta_e_g2")
    tree.Draw("eppi0.theta_g1_g2 * 180.0/TMath::Pi() >> h_theta_g1_g2")
    tree.Draw("eppi0.m2_miss >> h_m2_miss_pFD", "p.det==1")
    tree.Draw("eppi0.m2_miss >> h_m2_miss_pCD", "p.det==2")
    tree.Draw("eppi0.m2_epX >> h_m2_epX_pFD", "p.det==1")
    tree.Draw("eppi0.m2_epX >> h_m2_epX_pCD", "p.det==2")
    tree.Draw("eppi0.m2_epi0X >> h_m2_epi0X_pFD", "p.det==1")
    tree.Draw("eppi0.m2_epi0X >> h_m2_epi0X_pCD", "p.det==2")
    tree.Draw("eppi0.px_miss >> h_px_miss")
    tree.Draw("eppi0.py_miss >> h_py_miss")
    tree.Draw("eppi0.pz_miss >> h_pz_miss")
    tree.Draw("eppi0.E_miss >> h_E_miss_pFD", "p.det==1")
    tree.Draw("eppi0.E_miss >> h_E_miss_pCD", "p.det==2")
    tree.Draw("eppi0.pi0_deltaPhi * 180.0/TMath::Pi() >> h_deltaPhi")

    tree.Draw("dis.Q2 : dis.Xb >> h_Q2_vs_Xb")
    tree.Draw("dis.Q2 : dis.Xb >> h_Q2_vs_Xb_pFD", "p.det==1")
    tree.Draw("dis.Q2 : dis.Xb >> h_Q2_vs_Xb_pCD", "p.det==2")
    tree.Draw("eppi0.t : fmod((eppi0.trentoPhi * 180.0 / TMath::Pi()) + 360.0, 360.0) >> h_t_vs_phiT")
    tree.Draw("e.theta * 180.0/TMath::Pi() : e.p >> h_e_theta_vs_p")
    tree.Draw("e.phi * 180.0/TMath::Pi() : e.p >> h_e_phi_vs_p")
    tree.Draw("p.theta * 180.0/TMath::Pi() : p.p >> h_p_theta_vs_p")
    tree.Draw("p.phi * 180.0/TMath::Pi() : p.p >> h_p_phi_vs_p")
    tree.Draw("g.theta * 180.0/TMath::Pi() : g.p >> h_g_theta_vs_p")
    tree.Draw("g.phi * 180.0/TMath::Pi() : g.p >> h_g_phi_vs_p")
    tree.Draw("eppi0.m2_miss : p.theta * 180.0/TMath::Pi() >> h_m2_miss_vs_p_theta", "", "COLZ")
    tree.Draw("eppi0.m2_epX : p.theta * 180.0/TMath::Pi() >> h_m2_epX_vs_p_theta", "", "COLZ")
    tree.Draw("eppi0.m2_epi0X : p.theta * 180.0/TMath::Pi() >> h_m2_epi0X_vs_p_theta", "", "COLZ")
    tree.Draw("eppi0.m2_miss : eppi0.px_miss >> h_m2_miss_vs_px_miss", "", "COLZ")
    tree.Draw("eppi0.m2_miss : eppi0.py_miss >> h_m2_miss_vs_py_miss", "", "COLZ")
    tree.Draw("eppi0.m2_miss : eppi0.pz_miss >> h_m2_miss_vs_pz_miss", "", "COLZ")
    tree.Draw("eppi0.m2_miss : eppi0.E_miss >> h_m2_miss_vs_E_miss", "", "COLZ")

    outFile = ROOT.TFile.Open(outFilePath, "RECREATE")

    # ─── Fit Mgg histograms ─────────────────────
    def fit_mgg(hist):
        peak = hist.GetBinCenter(hist.GetMaximumBin())
        f = ROOT.TF1("gaus+pol1", "gaus(0)+pol1(3)", 0.098, 0.17)
        f.SetParameters(hist.GetMaximum(), peak, 0.01, 10, -1)
        hist.Fit(f, "RQ")
        # Make sure the fit gets drawn
        c = ROOT.TCanvas(f"c_{hist.GetName()}", "", 800, 600)
        hist.Draw()        # Draw histogram first
        f.Draw("SAME")     # Draw the fit on top
        c.Write()          # Save the canvas to the output file
        return f

    fit_mgg(histos1D["h_mgg_pFD"])
    fit_mgg(histos1D["h_mgg_pCD"])

    # ─── Save histograms & canvases ─────────────

    for h in histos1D.values():
        h.Write()

    for h in histos2D.values():
        h.Write()

    outFile.Close()
    print(f"All histograms and canvases saved to {outFilePath}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create EPPI0 exclusivity plots from ROOT TTree")
    parser.add_argument("input_file", type=str, help="Input ROOT file with TTree 'Events'")
    parser.add_argument("--output_file", type=str, help="Output ROOT file with histograms and fits")
    args = parser.parse_args()
    inputFile = args.input_file
    outputFile = outputFile = args.output_file or "EPPI0plots.root"
    EPPI0plots(inputFile, outputFile)
