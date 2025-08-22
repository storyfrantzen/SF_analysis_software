# type: ignore
import argparse
import ROOT
import numpy as np
import cProfile
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

def show_by_sector(particle, df, layer_name, y_branch, x_branch, sector_branch,
                        x_range=None, y_range=None, x_bins=100, y_bins=100, log_scale=False):
    """
    Plot 2D distributions by sector using ROOT histograms (multi-threaded).

    Parameters
    ----------
    particle : str
        Particle name for labeling
    df : ROOT.RDataFrame
        Input dataframe
    layer_name : str
        Layer identifier for labeling
    y_branch, x_branch : str
        Branch names to histogram
    sector_branch : str
        Branch identifying sector number
    x_range, y_range : tuple(float,float), optional
        Axis ranges
    x_bins, y_bins : int
        Histogram binning
    log_scale : bool
        If True, use log z-axis
    """

    ROOT.EnableImplicitMT()  # Make sure ROOT uses all cores

    canvas = ROOT.TCanvas(f"c_{particle}_{layer_name}", f"{particle} {layer_name}", 1200, 800)
    canvas.Divide(3, 2)

    hist_list = []

    for sec in range(1, 7):
        # Define histogram name/title
        hname = f"h_{particle}_{layer_name}_sec{sec}"
        title = f"{particle} {layer_name} Sector {sec}; {x_branch}; {y_branch}"

        # Fill ROOT 2D histogram
        hist = df.Filter(f"{sector_branch} == {sec}").Histo2D(
            (hname, title, x_bins, x_range[0], x_range[1], y_bins, y_range[0], y_range[1]),
            x_branch, y_branch
        )

        hist_list.append(hist)

        # Draw on canvas pad
        canvas.cd(sec)
        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetRightMargin(0.13)
        if log_scale:
            ROOT.gPad.SetLogz()
        hist.GetValue().Draw("COLZ")
        hist.GetValue().SetStats(False)

    canvas.Update()
    canvas.SaveAs(f"{particle}_{layer_name}_by_sector.png")

    return hist_list

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
        # show_by_sector("electron", df_ele, "ECAL", "SF", "rec.p", "rec.sector", x_range=(2, 7), y_range=(0.1, 0.35), x_bins=100, y_bins=100)
        # show_by_sector("photon",   df_pho, "PCAL", "rec.vPCAL", "rec.wPCAL", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        # show_by_sector("photon",   df_pho, "ECIN", "rec.vECIN", "rec.wECIN", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        # show_by_sector("photon",   df_pho, "ECOUT", "rec.vECOUT", "rec.wECOUT", "rec.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450))
        # show_by_sector("photon",   df_pho, "ECAL", "SF", "rec.p", "rec.sector", x_range=(0, 8), y_range=(0.1, 0.4), x_bins=100, y_bins=100)

        histos = show_by_sector(
        particle="electron",
        df=df_ele,
        layer_name="ECAL",
        y_branch="SF",
        x_branch="rec.p",
        sector_branch="rec.sector",
        x_range=(2, 7),
        y_range=(0.1, 0.35),
        x_bins=100,
        y_bins=100,
        )
        
    
    elif "e" in branch_list and "p" in branch_list and "g" in branch_list:
        # df = df.Define("eSF", "(e.E_PCAL + e.E_ECIN + e.E_ECOUT) / p")
        # df = df.Define("gSF", "(g.E_PCAL + g.E_ECIN + g.E_ECOUT) / p")
        # show_by_sector("electron", df, "PCAL", "e.vPCAL", "e.wPCAL", "e.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        # show_by_sector("electron", df, "ECIN", "e.vECIN", "e.wECIN", "e.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        # show_by_sector("electron", df, "ECOUT", "e.vECOUT", "e.wECOUT", "e.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        # show_by_sector("electron", df, "ECAL", "eSF", "e.p", "e.sector", x_range=(0, 8), y_range=(0, 0.4), x_bins=100, y_bins=100)
        # show_by_sector("photon",   df, "PCAL", "g.vPCAL", "g.wPCAL", "g.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        # show_by_sector("photon",   df, "ECIN", "g.vECIN", "g.wECIN", "g.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        # show_by_sector("photon",   df, "ECOUT", "g.vECOUT", "g.wECOUT", "g.sector", x_range=(0, 450), y_range=(0, 450), x_bins=450, y_bins=450)
        # show_by_sector("photon", df, "ECAL", "gSF", "g.p", "g.sector", x_range=(0, 8), y_range=(0, 0.4), x_bins=100, y_bins=100)
    
    else:
        print("Unable to load necessary particle branches! Check tree!")

if __name__ == "__main__":
    # Load BranchVars dictionary from install directory:
    ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")
    #cProfile.run("main()")
    main()