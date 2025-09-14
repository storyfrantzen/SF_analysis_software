# type: ignore
import argparse
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Macro computes the EPPI0 4-fold differential cross section from a skimmed EPPI0 root file

ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

def get_adaptive_edges(tree, varname, n_bins, min_val, max_val, fine_bins=500):
    """
    Adaptive binning based on equal counts per bin.
    
    Parameters:
        tree: ROOT TTree to draw from
        varname: variable name (string)
        n_bins: number of desired adaptive bins
        min_val, max_val: min and max values for histogram
        fine_bins: number of fine bins for initial histogram
    
    Returns:
        edges: list of bin edges (length = n_bins + 1)
    """
    # Create fine bin histogram
    hist_name = f"h_{varname}"
    hist = ROOT.TH1D(hist_name, hist_name, fine_bins, min_val, max_val)
    
    # Draw the variable into the histogram (without graphics)
    tree.Draw(f"{varname} >> {hist_name}", "", "goff")
    
    total_entries = hist.GetEntries()
    step = total_entries / n_bins
    
    edges = [hist.GetXaxis().GetXmin()]
    cumulative = 0
    target = step
    
    # Loop over fine bins to find edges where cumulative counts reach multiples of 'step'
    for i in range(1, hist.GetNbinsX() + 1):
        cumulative += hist.GetBinContent(i)
        if cumulative >= target:
            edge = hist.GetBinLowEdge(i + 1)  # upper edge of current bin
            if edge > edges[-1]:  # avoid duplicates
                edges.append(edge)
                target += step
    
    # Ensure the last edge is max_val
    if edges[-1] < hist.GetXaxis().GetXmax():
        edges.append(hist.GetXaxis().GetXmax())
    
    # Clean up histogram
    hist.Delete()
    
    return edges

### MAIN WORKFLOW: ###

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Compute XSection from input ROOT file.")
parser.add_argument("input_file", type=str, help="Path to the input EPPI0 ROOT file.")
parser.add_argument("--adaptive", action="store_true", help="Use adaptive binning scheme instead of fixed binning.")
parser.add_argument("--clas6", action="store_true", help="Use CLAS6 binning scheme.")

args = parser.parse_args()

# Get adaptive bin edges for each variable
if args.adaptive:
    Q2_edges = get_adaptive_edges(tree, "dis.Q2",   n_bins=4,  min_val=1.0,  max_val=6.5)
    Xb_edges = get_adaptive_edges(tree, "dis.Xb",   n_bins=4,  min_val=0.1,  max_val=0.7)
    t_edges  = get_adaptive_edges(tree, "eppi0.t",  n_bins=1,  min_val=0.0,  max_val=4)

    print("Adaptive Q2 edges:", Q2_edges)
    print("Adaptive Xb edges:", Xb_edges)
    print("Adaptive t edges:", t_edges)

if args.clas6:
    Q2_edges = np.array([1.0, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.6])
    Xb_edges = np.array([0.1, 0.15, 0.2, 0.25, 0.30, 0.38, 0.48, 0.58])
    t_edges = -1 * np.array([0.09, 0.15, 0.20, 0.30, 0.40, 0.60, 1.0, 1.5, 2.0])

# phi binning: uniform (for now). CLAS6 uses 20 bins.
phi_edges = np.linspace(0, 360, 21)  # final arg = bins + 1 

# Open ROOT file and access skimmed tree via dataframe:
file = ROOT.TFile(args.input_file, 'READ')
tree = file.Get('Events') 
df = ROOT.RDataFrame(tree)

# Snapshot needed columns into numpy arrays:
cols = ["dis.Q2", "dis.Xb", "eppi0.t", "eppi0.trentoPhi"]
df_numpy = df.AsNumpy(columns=cols)

Q2_vals = df_numpy["dis.Q2"]
Xb_vals = df_numpy["dis.Xb"]
t_vals = df_numpy["eppi0.t"]
phi_vals = (df_numpy["eppi0.trentoPhi"] + 2*np.pi) % (2*np.pi) * 180.0 / np.pi  # converted from radians to deg

# Helper function to find bin index
def find_bin(value, edges):
    idx = np.searchsorted(edges, value) - 1
    if idx < 0 or idx >= len(edges) - 1:
        return -1
    return idx
