import ROOT
ROOT.gROOT.SetBatch(True)

fData = ROOT.TFile.Open("RGKSKIM6.5_EPPI0plots.root")
fGEMC = ROOT.TFile.Open("PSEUDO6.5_EPPI0plots.root")
fOut = ROOT.TFile.Open("overlay_EPPI0plots_final.root", "RECREATE")

hist_names = [
    "h_helicity", "h_mgg", "h_mgg_pFD", "h_mgg_pCD",
    "h_Q2", "h_W", "h_t", "h_phiT",
    "h_theta_e_g1", "h_theta_e_g2", "h_theta_g1_g2",
    "h_m2_miss", "h_m2_miss_pFD", "h_m2_miss_pCD",
    "h_m2_epX", "h_m2_epX_pFD", "h_m2_epX_pCD",
    "h_m2_epi0X", "h_m2_epi0X_pFD", "h_m2_epi0X_pCD",
    "h_m_eggX", "h_m_eggX_pFD", "h_m_eggX_pCD",
    "h_px_miss", "h_py_miss", "h_pz_miss",
    "h_E_miss", "h_E_miss_pFD", "h_E_miss_pCD",
    "h_deltaPhi", "h_thetaX", "h_thetaX_pFD", "h_thetaX_pCD"
]

# Keep references to fits alive
all_funcs_to_keep = []

def extract_hist_and_funcs(obj):
    if not obj:
        return None, []
    h = None
    funcs = []
    if obj.InheritsFrom("TH1"):
        h = obj
        lof = obj.GetListOfFunctions()
        if lof:
            for f in lof:
                if f and f.InheritsFrom("TF1"):
                    funcs.append(f)
        return h, funcs
    if obj.InheritsFrom("TCanvas"):
        for prim in obj.GetListOfPrimitives():
            if prim.InheritsFrom("TH1"):
                h = prim
                lof = prim.GetListOfFunctions()
                if lof:
                    for f in lof:
                        if f and f.InheritsFrom("TF1"):
                            funcs.append(f)
                return h, funcs
    return None, []

def resample_hist_to_template(h_src, h_tmpl, new_name):
    nb = h_tmpl.GetNbinsX()
    xmin = h_tmpl.GetXaxis().GetXmin()
    xmax = h_tmpl.GetXaxis().GetXmax()
    h_new = ROOT.TH1D(new_name, h_src.GetTitle(), nb, xmin, xmax)
    for i in range(h_src.GetNbinsX() + 2):
        c = h_src.GetBinContent(i)
        if c:
            x = h_src.GetBinCenter(i) if 1 <= i <= h_src.GetNbinsX() else (xmin if i==0 else xmax)
            if x < xmin: x = xmin + 1e-9
            if x > xmax: x = xmax - 1e-9
            h_new.Fill(x, c)
    return h_new

for name in hist_names:
    hData, funcsData = extract_hist_and_funcs(fData.Get(name))
    hGEMC, funcsGEMC = extract_hist_and_funcs(fGEMC.Get(name))

    if not hData or not hGEMC:
        print(f"[skip] {name}: missing in one file")
        continue

    # Resample GEMC if needed
    same_binning = (hData.GetNbinsX() == hGEMC.GetNbinsX() and
                    abs(hData.GetXaxis().GetXmin() - hGEMC.GetXaxis().GetXmin()) < 1e-12 and
                    abs(hData.GetXaxis().GetXmax() - hGEMC.GetXaxis().GetXmax()) < 1e-12)
    if not same_binning:
        hGEMC = resample_hist_to_template(hGEMC, hData, hGEMC.GetName() + "_resampled")

    # Keep functions alive
    all_funcs_to_keep.extend(funcsData)
    all_funcs_to_keep.extend(funcsGEMC)

    # Styles
    hData.SetLineColor(ROOT.kBlack)
    hData.SetLineWidth(2)
    hData.SetStats(True)  # show one stats box for Data

    hGEMC.SetLineColor(ROOT.kRed)
    hGEMC.SetLineWidth(2)
    hGEMC.SetLineStyle(2)  # dashed histogram
    hGEMC.SetStats(True)   # show one stats box for GEMC

    ymax = max(hData.GetMaximum(), hGEMC.GetMaximum())
    hData.SetMaximum(ymax * 1.25)

    # Draw canvas
    c = ROOT.TCanvas(f"c_{name}", name, 900, 700)
    hData.Draw("hist")
    hGEMC.Draw("hist same")

    # Draw fits with solid lines
    for f in funcsData:
        if f:
            f.SetLineColor(ROOT.kBlack)
            f.SetLineWidth(2)
            f.SetLineStyle(1)  # solid
            f.Draw("same")
    for f in funcsGEMC:
        if f:
            f.SetLineColor(ROOT.kRed)
            f.SetLineWidth(2)
            f.SetLineStyle(1)  # solid
            f.Draw("same")

    # Legend
    leg = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.AddEntry(hData, "Data", "l")
    leg.AddEntry(hGEMC, "GEMC", "l")
    if funcsData: leg.AddEntry(funcsData[0], "Data fit", "l")
    if funcsGEMC: leg.AddEntry(funcsGEMC[0], "GEMC fit", "l")
    leg.Draw()

    # Write output
    fOut.cd()
    c.Write()
    hData.Write(name + "_data")
    hGEMC.Write(name + "_gemc")

fOut.Close()
print("overlay_EPPI0plots_final.root written.")
