# type: ignore
import argparse
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

def show_by_sector(particle, df, layer_name, y_branch, x_branch, sector_branch,
                   x_range=None, y_range=None, x_bins=100, y_bins=100, log_scale=False):

    fig, axes = plt.subplots(2, 3, figsize=(24, 9), constrained_layout=True)

    # Default bin ranges based on variable if not supplied
    if x_range is None:
        x_vals = df.AsNumpy([x_branch])[x_branch]
        x_range = (np.min(x_vals), np.max(x_vals))
    if y_range is None:
        y_vals = df.AsNumpy([y_branch])[y_branch]
        y_range = (np.min(y_vals), np.max(y_vals))

    bins_x = np.linspace(x_range[0], x_range[1], x_bins + 1)
    bins_y = np.linspace(y_range[0], y_range[1], y_bins + 1)

    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]

        # Filter for this sector
        df_sec = df.Filter(f"{sector_branch} == {sec}")
        arrays = df_sec.AsNumpy([x_branch, y_branch])
        x_vals = arrays[x_branch]
        y_vals = arrays[y_branch]

        # Build 2D histogram
        H, xedges, yedges = np.histogram2d(x_vals, y_vals, bins=[bins_x, bins_y])

        H_masked = np.ma.masked_where(H == 0, H)

        # Choose normalization
        if log_scale:
            norm = LogNorm(vmin=H_masked.min() if H_masked.min() > 0 else 1, vmax=H_masked.max())
        else:
            norm = Normalize(vmin=0, vmax=H_masked.max())

        # Create meshgrid from bin edges for pcolormesh
        X, Y = np.meshgrid(xedges, yedges)

        # Plot using pcolormesh
        pcm = ax.pcolormesh(X, Y, H_masked.T, cmap='viridis', norm=norm, shading='auto')

        ax.set_title(f"{layer_name} - Sector {sec}")
        ax.set_xlabel(f"{x_branch}")
        ax.set_ylabel(f"{y_branch}")
        fig.colorbar(pcm, ax=ax, label='Counts')

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
    #df = df.Range(int(5e5))

    # Get list of branch names
    branch_list = [b.GetName() for b in tree.GetListOfBranches()]

    if "rec" in branch_list:
        df = df.Define("SF", "(rec.E_PCAL + rec.E_ECIN + rec.E_ECOUT) / rec.p")
        df_ele = df.Filter("rec.pid == 11 && rec.charge == -1 && rec.sector > -1")
        df_pro = df.Filter("rec.pid == 2212 && rec.charge == 1")
        df_pho = df.Filter("rec.pid == 22 && rec.charge == 0")

        # show_by_sector("electron", df_ele, "PCAL", "rec.vPCAL", "rec.wPCAL", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        # show_by_sector("electron", df_ele, "ECIN", "rec.vECIN", "rec.wECIN", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        # show_by_sector("electron", df_ele, "ECOUT", "rec.vECOUT", "rec.wECOUT", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        show_by_sector("electron", df_ele, "ECAL", "SF", "rec.p", "rec.sector", x_range=(2, 7), y_range=(0.1, 0.35), x_bins=100, y_bins=100)
        # show_by_sector("photon",   df_pho, "PCAL", "rec.vPCAL", "rec.wPCAL", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        # show_by_sector("photon",   df_pho, "ECIN", "rec.vECIN", "rec.wECIN", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        # show_by_sector("photon",   df_pho, "ECOUT", "rec.vECOUT", "rec.wECOUT", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        show_by_sector("photon",   df_pho, "ECAL", "SF", "rec.p", "rec.sector", x_range=(0, 8), y_range=(0.1, 0.4), x_bins=100, y_bins=100)
    
    elif "e" in branch_list and "p" in branch_list and "g" in branch_list:
        df = df.Define("eSF", "(e.E_PCAL + e.E_ECIN + e.E_ECOUT) / p")
        df = df.Define("gSF", "(g.E_PCAL + g.E_ECIN + g.E_ECOUT) / p")
        show_by_sector("electron", df, "PCAL", "e.vPCAL", "e.wPCAL", "e.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        show_by_sector("electron", df, "ECIN", "e.vECIN", "e.wECIN", "e.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        show_by_sector("electron", df, "ECOUT", "e.vECOUT", "e.wECOUT", "e.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        show_by_sector("electron", df, "ECAL", "eSF", "e.p", "e.sector", x_range=(0, 8), y_range=(0, 0.4), x_bins=100, y_bins=100, log_scale=True)
        show_by_sector("photon",   df, "PCAL", "g.vPCAL", "g.wPCAL", "g.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        show_by_sector("photon",   df, "ECIN", "g.vECIN", "g.wECIN", "g.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        show_by_sector("photon",   df, "ECOUT", "g.vECOUT", "g.wECOUT", "g.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        show_by_sector("photon", df, "ECAL", "gSF", "g.p", "g.sector", x_range=(0, 8), y_range=(0, 0.4), x_bins=100, y_bins=100)
    
    else:
        print("Unable to load necessary particle branches! Check tree!")

if __name__ == "__main__":
    # Load BranchVars dictionary from install directory:
    ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")
    main()