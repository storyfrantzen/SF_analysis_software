# chi2pid_vs_p.py
# type: ignore
import ROOT

# Open ROOT file
file = ROOT.TFile.Open("output/6.5inclusiveRecon50_noCuts.root")
tree = file.Get("Events")

# Create canvas
c = ROOT.TCanvas("c", "proton beta vs p", 800, 600)

# Define 2D histogram
hist = ROOT.TH2D("h_p_beta_vs_p_p", ";p [GeV]; #beta", 
                 100, 0, 6.5,    # x-axis: p_p
                 100, 0, 1)     # y-axis: chi2pid (adjust range if needed)

# Fill histogram using TTree::Draw
tree.Draw("p_beta : p_p >> h_p_beta_vs_p_p", "colz")

# Draw and save
hist.Draw("colz")
c.SaveAs("p_beta_vs_p_p.png")
