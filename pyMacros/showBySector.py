# type: ignore
import argparse
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def show_by_sector(particle, df, layer_name, y_branch, x_branch, sector_branch):
    fig, axes = plt.subplots(2, 3, figsize=(18, 9), constrained_layout=True)
    bins_x = np.linspace(0, 450, 451)
    bins_y = np.linspace(0, 450, 451)

    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]

        # Filter for this sector
        df_sec = df.Filter(f"{sector_branch} == {sec}")

        # Directly get NumPy arrays from RDataFrame
        arrays = df_sec.AsNumpy([x_branch, y_branch])
        x_vals = arrays[x_branch]
        y_vals = arrays[y_branch]

        # Build 2D histogram
        H, xedges, yedges = np.histogram2d(x_vals, y_vals, bins=[bins_x, bins_y])

        # Plot
        im = ax.imshow(
            H.T, origin='lower',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            aspect='auto', interpolation='nearest', cmap='viridis',
            norm=LogNorm(vmin=1, vmax=H.max() if H.max() > 0 else 1)
        )

        ax.set_title(f"{layer_name} - Sector {sec}")
        ax.set_xlabel(f"{x_branch} [cm]")
        ax.set_ylabel(f"{y_branch} [cm]")
        fig.colorbar(im, ax=ax, orientation='vertical', label='Counts')

    plt.suptitle(f"2D Histograms of {layer_name}: {y_branch} vs {x_branch} by Sector")
    plt.savefig(f"{particle}_{layer_name}_bySector.png")
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Create acceptance plots from ROOT TTree using RDataFrame and matplotlib")
    parser.add_argument("input_file", type=str, help="Input ROOT file with TTree 'Events'")
    args = parser.parse_args()

    infile = ROOT.TFile.Open(args.input_file)
    if not infile or infile.IsZombie():
        print(f"Error: Could not open input file {args.input_file}")
        return

    tree = infile.Get("Events")
    if not tree:
        print(f"Error: 'Events' tree not found in {args.input_file}")
        return

    df = ROOT.RDataFrame(tree)

     # Get list of branch names
    branch_list = [b.GetName() for b in tree.GetListOfBranches()]

    if "rec" in branch_list:
        df_ele = df.Filter("rec.pid == 11 && rec.charge == -1 && rec.sector > -1")
        df_pro = df.Filter("rec.pid == 2212 && rec.charge == 1")
        df_pho = df.Filter("rec.pid == 22 && rec.charge == 0")

        show_by_sector("electron", df_ele, "PCAL", "rec.vPCAL", "rec.wPCAL", "rec.sector")
        show_by_sector("electron", df_ele, "ECIN", "rec.vECIN", "rec.wECIN", "rec.sector")
        show_by_sector("electron", df_ele, "ECOUT", "rec.vECOUT", "rec.wECOUT", "rec.sector")
        show_by_sector("photon",   df_pho, "PCAL", "rec.vPCAL", "rec.wPCAL", "rec.sector")
        show_by_sector("photon",   df_pho, "ECIN", "rec.vECIN", "rec.wECIN", "rec.sector")
        show_by_sector("photon",   df_pho, "ECOUT", "rec.vECOUT", "rec.wECOUT", "rec.sector")
    
    elif "e" in branch_list and "p" in branch_list and "g" in branch_list:

        show_by_sector("electron", df, "PCAL", "e.vPCAL", "e.wPCAL", "e.sector")
        show_by_sector("electron", df, "ECIN", "e.vECIN", "e.wECIN", "e.sector")
        show_by_sector("electron", df, "ECOUT", "e.vECOUT", "e.wECOUT", "e.sector")
        show_by_sector("photon",   df, "PCAL", "g.vPCAL", "g.wPCAL", "g.sector")
        show_by_sector("photon",   df, "ECIN", "g.vECIN", "g.wECIN", "g.sector")
        show_by_sector("photon",   df, "ECOUT", "g.vECOUT", "g.wECOUT", "g.sector")
    
    else:
        print("Unable to load necessary particle branches! Check tree!")

if __name__ == "__main__":
    # Load BranchVars dictionary from install directory:
    ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")
    main()