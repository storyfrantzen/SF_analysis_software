
#type: ignore
import ROOT
import argparse
import numpy as np

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Compute proton energy loss corrections from input MC ROOT file.")
parser.add_argument("input_file", type=str, help="Path to the input ROOT file")
parser.add_argument("--isCD", action="store_true", help="Proton detected in CD")
args = parser.parse_args()

# Open MC ROOT file
file = ROOT.TFile(args.input_file, 'READ')

tree = file.Get("Events")  # or whatever your TTree is named

# Define theta binning
theta_edges = np.linspace(35, 70, 11) if args.isCD else np.linspace(12, 45, 25)

# Create output canvas
n_bins = len(theta_edges) - 1
n_rows = 2 if args.CD else 4
n_cols = n_bins / n_rows
# n_rows = (n_bins + n_cols - 1) // n_cols  # Ceiling division
canvas = ROOT.TCanvas("c_p_deltap_vs_p_p", "Proton Delta p vs p in theta bins", 2000, 350 * n_rows)
canvas.Divide(n_cols, n_rows)

# Create output file for histograms
output_file = ROOT.TFile("proton_deltap_thetaBins.root", "RECREATE")

# Store histograms in a list
hist_list = []

# Loop over theta bins
for i in range(n_bins):
    tmin = theta_edges[i]
    tmax = theta_edges[i + 1]
    
    # Define histogram name
    hname = f"h_p_deltap_p_theta_{int(tmin)}_{int(tmax)}"

    if args.isCD:
        y_range = (-0.2, 0.2)
    else:
        y_range = (-0.065, 0.2) if tmin > 35 else (-0.065, 0.065)
    hist = ROOT.TH2F(hname, f"Proton #Delta p vs p, {tmin:.1f} < #theta < {tmax:.1f}; " + "p_{rec} [GeV]; #Delta p [GeV]",
                100, 0, 5, 100, y_range[0], y_range[1])

    # Build draw command and cut
    draw_cmd = "(p_pgen - p_p) : p_p >> " + hname
    degtoRad = np.pi / 180

    if args.isCD:
        cut = f"p_theta > {tmin * degtoRad} && p_theta < {tmax * degtoRad} && pid == 2212 && detPro == 2"
    else:
        cut = f"p_theta > {tmin * degtoRad} && p_theta < {tmax * degtoRad} && pid == 2212 && detPro == 1"

    # Fill histogram
    tree.Draw(draw_cmd, cut, "goff")

    # Write to file
    hist.Write()
    hist_list.append(hist)

    # Draw on canvas
    canvas.cd(i + 1)
    pad = ROOT.gPad
    pad.SetLeftMargin(0.15)   # Increase to give room for y-axis label
    pad.SetLogz()       # Set log scale on z-axis
    hist.Draw("COLZ")

    # Create profile from the 2D histogram
    prof = hist.ProfileX()
    prof.SetLineColor(ROOT.kBlack)
    prof.SetMarkerColor(ROOT.kBlack)
    prof.SetMarkerStyle(10)

    nbins = prof.GetNbinsX()
    xmin = None
    xmax = None

    for b in range(1, nbins + 1):
        if prof.GetBinEntries(b) > 0:
            x = prof.GetBinCenter(b)
            if xmin is None:
                xmin = x  # first non-empty bin
            xmax = x      # keep updating as we go right

    if xmin is not None and xmax is not None and xmax > xmin:
        if args.isCD:
            fit_func = ROOT.TF1(f"f_fit_{i}", "[0] + [1] * x + [2] * (x*x)", xmin, xmax)
        else:
            fit_func = ROOT.TF1(f"f_fit_{i}", "[0] + [1] / x + [2] / (x*x)", xmin, xmax)
        prof.Fit(fit_func, "R")  # 'R' ensures fit only within [xmin, xmax]
        fit_func.SetLineColor(ROOT.kRed)
        fit_func.SetLineWidth(1) 

    # 4. Now re-draw the 2D hist axes and colormap (without clearing pad)
    hist.Draw("COLZ")  # Re-assert 2D histogram

    # 5. Overlay profile and fit
    prof.Draw("SAME")
    fit_func.Draw("SAME")

    # Draw the profile and the fit over the 2D histogram

# Save canvas
canvas.SaveAs("proton_deltap_vs_p_thetaBins.png")
output_file.Close()
print("Histograms and canvas saved.")
