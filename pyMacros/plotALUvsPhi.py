# type: ignore
import ROOT
import numpy as np
import matplotlib.pyplot as plt

n_bins = 12
phi_min = 0
phi_max = 360
P = 0.8

# Open the ROOT file
file = ROOT.TFile('output/6.5eppi0Skim18_allfid.root', 'READ')

# Access the tree inside the file
tree = file.Get('Events') 

rdf = ROOT.RDataFrame(tree)
rdf = rdf.Define("phi_deg", "fmod(trentoPhi + 2*TMath::Pi(), 2*TMath::Pi()) * 180.0 / TMath::Pi()")
rdf_pos = rdf.Filter("helicity == 1")
rdf_neg = rdf.Filter("helicity == -1")

canvas = ROOT.TCanvas("c1", "Phi Distributions", 800, 600)

hist_phi_pos = rdf_pos.Histo1D(("hpos_phi", "phi (deg)", n_bins, phi_min, phi_max), "phi_deg")
hist_phi_pos.SetLineColor(ROOT.kRed)
hist_phi_pos.SetLineWidth(2)
hist_phi_pos.SetTitle("Phi Distribution for Positive Helicity")
hist_phi_pos.Draw()

hist_phi_neg = rdf_neg.Histo1D(("hneg_phi", "phi (deg)", n_bins, phi_min, phi_max), "phi_deg")
hist_phi_neg.SetLineColor(ROOT.kBlue)
hist_phi_neg.SetLineWidth(2)
hist_phi_neg.Draw("Same")

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(hist_phi_pos.GetPtr(), "Positive Helicity", "l")
legend.AddEntry(hist_phi_neg.GetPtr(), "Negative Helicity", "l")
legend.Draw()

hist_sum = hist_phi_pos.Clone("hist_sum") # Clone the positive helicity histogram to preserve it
hist_sum.SetTitle("Phi Sum (Positive + Negative)")
hist_sum.SetLineColor(ROOT.kBlue)
hist_sum.SetLineWidth(2)

hist_sum.Add(hist_phi_neg.GetValue())

hist_diff = hist_phi_pos.Clone("hist_diff")  
hist_diff.SetTitle("Phi Difference (Positive - Negative)")
hist_diff.SetLineColor(ROOT.kBlack)
hist_diff.SetLineWidth(2)

hist_diff.Add(hist_phi_neg.GetValue(), -1) # subtract negative helicity hist, dereferencing the pointer to the negative hist
hist_diff.Divide(hist_sum)

hist_Asym = hist_diff.Clone("hist_Asym")
hist_Asym.Scale(1.0 / P)

for i in range(1, hist_Asym.GetNbinsX() + 1):
    N_plus = hist_phi_pos.GetBinContent(i)
    N_minus = hist_phi_neg.GetBinContent(i)
    N_total = hist_sum.GetBinContent(i)

    if N_total > 0:
        # Compute the asymmetry bin error using the provided formula
        asymmetry_error = (2 / P) * np.sqrt((N_plus * N_minus) / (N_total ** 3))
        hist_Asym.SetBinError(i, asymmetry_error)
    else:
        # Set error to zero if no events in the bin
        hist_Asym.SetBinError(i, 0)

# Define a sine function for the fit
sin_func = ROOT.TF1("sin_func", "[0] * sin((x + [1]) * TMath::DegToRad())", phi_min, phi_max)
sin_func.SetParameters(0, 1, 0)  # Initial guess for parameters: offset, amplitude, phase shift
hist_Asym.Fit(sin_func, "ER")  # Fit with the sine function

ROOT.gStyle.SetErrorX(0.0001)
ROOT.gStyle.SetOptStat(0)

hist_Asym.SetMarkerStyle(20)
hist_Asym.Draw("E1")
sin_func.Draw("same") 

hist_Asym.SetTitle("$\\frac{N_+ - N_-}{N_+ + N_-}$ vs $\\phi$")
hist_Asym.GetXaxis().SetTitle("$\\Phi$")
hist_Asym.GetYaxis().SetTitle("Counts")

canvas.Update()
canvas.SaveAs("BSA vs Phi.png")