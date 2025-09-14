#type: ignore
import ROOT
import argparse
import sys

def fit_histogram(hist, func_str, fit_range, init_params=None, outFile=None):
    """
    Fits a histogram with the specified function, draws it with ±3σ lines,
    and saves the canvas into the provided ROOT file.

    Parameters
    ----------
    hist : ROOT.TH1
        Histogram to fit.
    func_str : str
        Fit function string, e.g. "gaus(0)+pol1(3)".
    fit_range : tuple
        (min, max) range for the fit.
    init_params : list, optional
        Initial parameters for the fit.
    outFile : ROOT.TFile, optional
        ROOT file to save the canvas to.
    """
    # Determine approximate peak
    peak = hist.GetBinCenter(hist.GetMaximumBin())
    
    # Create fit function
    f = ROOT.TF1(f"f_{hist.GetName()}", func_str, *fit_range)

    # Set initial parameters
    if init_params:
        for i, val in enumerate(init_params):
            f.SetParameter(i, val)
    else:
        # Default example for gaus(0)+pol1(3) (6 params)
        f.SetParameters(hist.GetMaximum(), peak, 0.01, 0, 0, 0)

    # Perform fit quietly
    hist.Fit(f, "RQ")

    # Create canvas
    c = ROOT.TCanvas(f"c_{hist.GetName()}", f"Fit {hist.GetName()}", 800, 600)
    hist.Draw()      # Draw histogram
    f.Draw("SAME")   # Draw fit function on top

    # ±3σ lines
    mean = f.GetParameter(1)
    sigma = f.GetParameter(2)
    lines = []  # keep references
    for shift in [-3, 3]:
        line = ROOT.TLine(mean + shift*sigma, 0, mean + shift*sigma, hist.GetMaximum())
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)  # dashed
        line.SetLineWidth(2)
        line.Draw("SAME")
        lines.append(line) # keep reference!
    
    # Mean line (blue solid)
    mean_line = ROOT.TLine(mean, 0, mean, hist.GetMaximum())
    mean_line.SetLineColor(ROOT.kBlue)
    mean_line.SetLineStyle(2)
    mean_line.SetLineWidth(1)
    mean_line.Draw("SAME")
    lines.append(mean_line)

    c.Update()

    # Save canvas to ROOT file if provided
    if outFile:
        outFile.cd()
        c.Write()

    return f


# ─── Specialized Fits ─────────────────────
def fit_mgg(hist, outFile=None):
    return fit_histogram(hist, "gaus(0)+pol1(3)", (0.098, 0.17),
                            init_params=[hist.GetMaximum(), 0.135, 0.01, 10, -1], outFile=outFile)

def fit_m2_epX(hist, outFile=None):
    return fit_histogram(hist, "gaus(0)+pol1(3)", (-0.2, 0.2), outFile=outFile)

def fit_m2_epi0X(hist, outFile=None):
    return fit_histogram(hist, "gaus(0)+pol1(3)", (0.5, 1.5), outFile=outFile)

def fit_m_epi0X(hist, outFile=None):
    return fit_histogram(hist, "gaus(0)+pol1(3)", (0.8, 1.1), outFile=outFile)

def fit_E_miss(hist, outFile=None):
    return fit_histogram(hist, "gaus(0)+pol1(3)", (-0.2, 0.4), outFile=outFile)


def EPPI0plots(inputFilePath="input.root", outFilePath="EPPI0plots.root"):
    ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

    f = ROOT.TFile.Open(inputFilePath)
    if not f or f.IsZombie():
        print(f"Could not open file: {inputFilePath}")
        return

    tree = f.Get("Events")
    if not tree:
        print("Could not find tree 'Events' in file.")
        return

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
        "h_m_epi0X_pFD": ROOT.TH1D("h_m_epi0X_pFD", "pFD M_{e#piX}; M_{e#piX} [GeV]; Counts", 200, -1, 4),
        "h_m_epi0X_pCD": ROOT.TH1D("h_m_epi0X_pCD", "pCD M_{e#piX}; M_{e#piX} [GeV]; Counts", 200, -1, 4),
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
    #tree.Draw("sqrt(eppi0.m2_epi0X) >> h_m_epi0X_pFD", "p.det==1")
    #tree.Draw("sqrt(eppi0.m2_epi0X) >> h_m_epi0X_pCD", "p.det==2")

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


    # # ─── Run Fits for pFD and pCD ─────────────────────
    fits = {}
    results = {}

    # Mgg
    fits["mgg_pFD"]    = fit_mgg(histos1D["h_mgg_pFD"], outFile=outFile)
    fits["mgg_pCD"]    = fit_mgg(histos1D["h_mgg_pCD"], outFile=outFile)

    # MM²(epX)
    fits["m2_epX_pFD"] = fit_m2_epX(histos1D["h_m2_epX_pFD"], outFile=outFile)
    fits["m2_epX_pCD"] = fit_m2_epX(histos1D["h_m2_epX_pCD"], outFile=outFile)

    # MM²(epi0X)
    fits["m2_epi0X_pFD"] = fit_m2_epi0X(histos1D["h_m2_epi0X_pFD"], outFile=outFile)
    fits["m2_epi0X_pCD"] = fit_m2_epi0X(histos1D["h_m2_epi0X_pCD"], outFile=outFile)
    

    # # MX(eγγ)
    # #fits["m_epi0X_pFD"] = fit_m_epi0X(histos1D["h_m_epi0X_pFD"])
    # #fits["m_epi0X_pCD"] = fit_m_epi0X(histos1D["h_m_epi0X_pCD"])

    # Emiss
    fits["E_miss_pFD"]  = fit_E_miss(histos1D["h_E_miss_pFD"], outFile=outFile)
    fits["E_miss_pCD"]  = fit_E_miss(histos1D["h_E_miss_pCD"], outFile=outFile)

    # ─── Extract Mean and Sigma for Later Cuts ─────────────────────
    # for key, f in fits.items():
    #     if f:  # make sure the fit is valid
    #         mean  = f.GetParameter(1)
    #         sigma = f.GetParameter(2)
    #         results[key] = (mean, sigma)
    #     else:
    #         print(f"Warning: Fit for {key} failed.")
    #         results[key] = (None, None)

    # print("Exclusivity fit results:")
    # for key, (mean, sigma) in results.items():
    #     print(f"{key}: mean = {mean:.4f}, sigma = {sigma:.4f}")

    # ─── Save histograms & canvases ─────────────

    for h in histos1D.values(): 
        if h is not None:
            h.Write()

    for h in histos2D.values(): 
        if h is not None:
            h.Write()

    outFile.Close()
    print(f"All histograms and canvases saved to {outFilePath}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create EPPI0 exclusivity plots from ROOT TTree")
    parser.add_argument("input_file", type=str, help="Input ROOT file with TTree 'Events'")
    parser.add_argument("--output_file", type=str, help="Output ROOT file with histograms and fits")
    args = parser.parse_args()
    inputFile = args.input_file
    outputFile = args.output_file or "EPPI0plots.root"
    EPPI0plots(inputFile, outputFile)
