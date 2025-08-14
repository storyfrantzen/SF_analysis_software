# type: ignore
import os
import argparse
import ROOT

ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

ROOT.gInterpreter.Declare("""
double wrap360(double deg) {
  return (deg < 0) ? deg + 360.0 : deg;
}
""")
def drawElectronPlots(df, branchPrefix, leaf_list):
    # Mandatory variables
    p_var = f"{branchPrefix}.p"
    beta_var = f"{branchPrefix}.beta"
    theta_var = f"{branchPrefix}_theta_deg"
    phi_var = f"{branchPrefix}_phi_deg"
    vz_var = f"{branchPrefix}.vz"
    sector_var = f"{branchPrefix}.sector"
    det_var = f"{branchPrefix}.det"

    # Optional DC coordinates
    xDC1_var = f"{branchPrefix}.xDC1"
    yDC1_var = f"{branchPrefix}.yDC1"
    xDC2_var = f"{branchPrefix}.xDC2"
    yDC2_var = f"{branchPrefix}.yDC2"
    xDC3_var = f"{branchPrefix}.xDC3"
    yDC3_var = f"{branchPrefix}.yDC3"

    e_hists = []

    # Always-make histograms
    e_hists.append(df.Histo1D(("e_p", "Electron Momentum; p_{e}; Counts", 300, 0, 6.5), p_var))
    e_hists.append(df.Histo1D(("e_beta", "Electron Beta; beta_{e}; Counts", 100, 0, 2), beta_var))
    e_hists.append(df.Histo1D(("e_theta", "Electron Theta; #theta_{e} [deg]; Counts", 120, 0, 120), theta_var))
    e_hists.append(df.Histo1D(("e_phi", "Electron Phi; #phi_{e} [deg]; Counts", 360, 0, 360), phi_var))
    e_hists.append(df.Histo1D(("e_vz", "Electron v_{z}; #v_{z} [cm]; Counts", 20, -10, 10), vz_var))
    e_hists.append(df.Histo1D(("e_sector", "Electron Sector; Sector; Counts", 6, 1, 7), sector_var))
    e_hists.append(df.Histo1D(("e_det", "Electron Detector; Detector; Counts", 3, 0, 3), det_var))
    e_hists.append(df.Histo2D(("e_beta_vs_p", "Electron #beta vs p; p [GeV]; #beta", 500, 0, 10, 500, 0, 2), p_var, beta_var))

    # Only make DC histograms if all required branches exist
    if all(v in leaf_list for v in ["xDC1", "yDC1"]):
        e_hists.append(df.Histo2D(("e_yDC1_vs_xDC1", "Electron DC Region 1 y vs x; x [cm]; y [cm]",
                                 300, -150, 150, 300, -150, 150), xDC1_var, yDC1_var))

    if all(v in leaf_list for v in ["xDC2", "yDC2"]):
        e_hists.append(df.Histo2D(("e_yDC2_vs_xDC2", "Electron DC Region 2 y vs x; x [cm]; y [cm]",
                                 400, -200, 200, 400, -200, 200), xDC2_var, yDC2_var))

    if all(v in leaf_list for v in ["xDC3", "yDC3"]):
        e_hists.append(df.Histo2D(("e_yDC3_vs_xDC3", "Electron DC Region 3 y vs x; x [cm]; y [cm]",
                                 400, -300, 300, 400, -300, 300), xDC3_var, yDC3_var))

    return e_hists

def drawProtonPlots(df, branchPrefix, leaf_list):
    p_var = f"{branchPrefix}.p"
    beta_var = f"{branchPrefix}.beta"
    theta_var = f"{branchPrefix}_theta_deg"  
    phi_var = f"{branchPrefix}_phi_deg"      
    vz_var = f"{branchPrefix}.vz"
    sector_var = f"{branchPrefix}.sector"
    det_var = f"{branchPrefix}.det"

    # Optional DC coordinates
    xDC1_var = f"{branchPrefix}.xDC1"
    yDC1_var = f"{branchPrefix}.yDC1"
    xDC2_var = f"{branchPrefix}.xDC2"
    yDC2_var = f"{branchPrefix}.yDC2"
    xDC3_var = f"{branchPrefix}.xDC3"
    yDC3_var = f"{branchPrefix}.yDC3"

    # Optional CVT coordinates
    edge_cvt1_var = f"{branchPrefix}.edge_cvt1"
    edge_cvt3_var = f"{branchPrefix}.edge_cvt3"
    edge_cvt5_var = f"{branchPrefix}.edge_cvt5"
    edge_cvt7_var = f"{branchPrefix}.edge_cvt7"
    edge_cvt12_var = f"{branchPrefix}.edge_cvt12"
    theta_cvt_var = f"{branchPrefix}_theta_cvt_deg"
    phi_cvt_var = f"{branchPrefix}_phi_cvt_deg"

    p_hists = []

    p_hists.append(df.Histo1D(("p_p", "Proton Momentum; p_{p}; Counts", 300, 0, 6.5), p_var))
    p_hists.append(df.Histo1D(("p_beta", "Proton Beta; beta_{p}; Counts", 50, 0, 2), beta_var))
    p_hists.append(df.Histo1D(("p_theta", "Proton Theta; #theta_{p} [deg]; Counts", 120, 0, 120), theta_var))
    p_hists.append(df.Histo1D(("p_phi", "Proton Phi; #phi_{p} [deg]; Counts", 360, 0, 360), phi_var))
    p_hists.append(df.Histo1D(("p_vz", "Proton v_{z}; #v_{z} [cm]; Counts", 20, -10, 10), vz_var))
    p_hists.append(df.Histo1D(("p_sector", "Proton Sector; Sector; Counts", 6, 1, 7), sector_var))
    p_hists.append(df.Histo1D(("p_det", "Proton Detector; Detector; Counts", 3, 0, 3), det_var))
    p_hists.append(df.Histo2D(("p_beta_vs_p", "Proton #beta vs p; p [GeV]; #beta", 500, 0, 10, 500, 0, 2), p_var, beta_var))

    # Only make DC histograms if all required branches exist
    if all(v in leaf_list for v in ["xDC1", "yDC1"]):
        p_hists.append(df.Histo2D(("p_yDC1_vs_xDC1", "Proton DC Region 1 y vs x; x [cm]; y [cm]",
                                 300, -150, 150, 300, -150, 150), xDC1_var, yDC1_var))

    if all(v in leaf_list for v in ["xDC2", "yDC2"]):
        p_hists.append(df.Histo2D(("p_yDC2_vs_xDC2", "Proton DC Region 2 y vs x; x [cm]; y [cm]",
                                 400, -200, 200, 400, -200, 200), xDC2_var, yDC2_var))

    if all(v in leaf_list for v in ["xDC3", "yDC3"]):
        p_hists.append(df.Histo2D(("p_yDC3_vs_xDC3", "Proton DC Region 3 y vs x; x [cm]; y [cm]",
                                 400, -300, 300, 400, -300, 300), xDC3_var, yDC3_var))
        
    
    if all(v in leaf_list for v in ["phi_cvt", "theta_cvt"]):
        p_hists.append(df.Histo2D(("p_theta_cvt_vs_phi_cvt", "Proton CVT #theta vs #phi; #phi [deg]; #theta [deg]",
                                   360, 0, 360, 140, 0, 140), phi_cvt_var, theta_cvt_var))

    return p_hists

def drawPhotonPlots(df, branchPrefix, leaf_list):
    p_var = f"{branchPrefix}.p"
    beta_var = f"{branchPrefix}.beta"
    theta_var = f"{branchPrefix}_theta_deg"  
    phi_var = f"{branchPrefix}_phi_deg"      
    vz_var = f"{branchPrefix}.vz"
    sector_var = f"{branchPrefix}.sector"
    det_var = f"{branchPrefix}.det"

    g_hists = []

    g_hists.append(df.Histo1D(("g_p", "Photon Momentum; p_{#gamma}; Counts", 300, 0, 6.5), p_var))
    g_hists.append(df.Histo1D(("g_beta", "Photon Beta; beta_{#gamma}; Counts", 50, 0, 2), beta_var))
    g_hists.append(df.Histo1D(("g_theta", "Photon Theta; #theta_{#gamma} [deg]; Counts", 120, 0, 120), theta_var))
    g_hists.append(df.Histo1D(("g_phi", "Photon Phi; #phi_{#gamma} [deg]; Counts", 360, 0, 360), phi_var))
    g_hists.append(df.Histo1D(("g_vz", "Photon v_{z}; #v_{z} [cm]; Counts", 20, -10, 10), vz_var))
    g_hists.append(df.Histo1D(("g_sector", "Photon Sector; Sector; Counts", 6, 1, 7), sector_var))
    g_hists.append(df.Histo1D(("g_det", "Photon Detector; Detector; Counts", 3, 0, 3), det_var))
    g_hists.append(df.Histo2D(("g_beta_vs_p", "Photon #beta vs p; p [GeV]; #beta", 500, 0, 10, 500, 0, 2), p_var, beta_var))

    return g_hists


def main():
    parser = argparse.ArgumentParser(description="Create acceptance plots from ROOT TTree using RDataFrame")
    parser.add_argument("input_file", type=str, help="Input ROOT file with TTree 'Events'")
    parser.add_argument("--max_events", type=int, help="Maximum number of events to visualize")
    args = parser.parse_args()

    # Open input file and get tree
    infile = ROOT.TFile.Open(args.input_file)
    if not infile or infile.IsZombie():
        print(f"Error: Could not open input file {args.input_file}")
        return
    tree = infile.Get("Events")
    if not tree:
        print(f"Error: 'Events' tree not found in {args.input_file}")
        return

    df = ROOT.RDataFrame(tree)
    if args.max_events is not None:
        df = df.Range(args.max_events)
    
    leaf_list = [leaf.GetName() for leaf in tree.GetListOfLeaves()]
    #print(leaf_list)

    all_hists = []

    if tree.GetBranch("rec"):
        print("Input file is a rec file. Proceed with branch rec.")
        # Define new columns for degrees (use valid C++ expressions)
        df = df.Define("rec_theta_deg", "wrap360(rec.theta * 180.0 / 3.141592653589793)")
        df = df.Define("rec_phi_deg", "wrap360(rec.phi * 180.0 / 3.141592653589793)")
        if all(v in leaf_list for v in ["phi_cvt", "theta_cvt"]):
            df = df.Define("rec_theta_cvt_deg", "wrap360(rec.theta_cvt * 180.0 / 3.141592653589793)")
            df = df.Define("rec_phi_cvt_deg", "wrap360(rec.phi_cvt * 180.0 / 3.141592653589793)")
        # Filters for particle types by pid
        df_ele = df.Filter("rec.pid == 11 && rec.charge == -1 && rec.sector > -1")
        df_pro = df.Filter("rec.pid == 2212 && rec.charge == 1")
        df_pho = df.Filter("rec.pid == 22 && rec.charge == 0")

        h_pid = df.Histo1D(("pid", "Reconstructed particle PID; PID; Counts", 500, -3000, 3000), "rec.pid")
        h_beta_vs_p = df.Histo2D(("beta_vs_p", "Reconstructed #beta vs p; p [GeV]; #beta", 500, 0, 10, 500, 0, 2), "rec.p", "rec.beta")

        all_hists.append(h_pid)
        all_hists.append(h_beta_vs_p)

        e_hists = drawElectronPlots(df_ele, "rec", leaf_list)
        p_hists = drawProtonPlots(df_pro, "rec", leaf_list)
        g_hists = drawPhotonPlots(df_pho, "rec", leaf_list)
    
    elif tree.GetBranch("e") and tree.GetBranch("p") and tree.GetBranch("g"):
        print("Input file is an EPPI0 file. Proceed with branches e, p, and g.")
        df = df.Define("e_theta_deg", "wrap360(e.theta * 180.0 / 3.141592653589793)")
        df = df.Define("e_phi_deg", "wrap360(e.phi * 180.0 / 3.141592653589793)")
        df = df.Define("p_theta_deg", "wrap360(p.theta * 180.0 / 3.141592653589793)")
        df = df.Define("p_phi_deg", "wrap360(p.phi * 180.0 / 3.141592653589793)")
        if all(v in leaf_list for v in ["phi_cvt", "theta_cvt"]):
            print("setting p_theta_cvt_deg && p_phi_cvt")
            df = df.Define("p_theta_cvt_deg", "wrap360(p.theta_cvt * 180.0 / 3.141592653589793)")
            df = df.Define("p_phi_cvt_deg", "wrap360(p.phi_cvt * 180.0 / 3.141592653589793)")
        df = df.Define("g_theta_deg", "wrap360(g.theta * 180.0 / 3.141592653589793)")
        df = df.Define("g_phi_deg", "wrap360(g.phi * 180.0 / 3.141592653589793)")

        e_hists = drawElectronPlots(df, "e", leaf_list)
        p_hists = drawProtonPlots(df, "p", leaf_list)
        g_hists = drawPhotonPlots(df, "g", leaf_list)
    
    else:
        print("Tree does not have expected branches 'rec' or ('e','p','g'). Exiting.")
        return
    
    for hist in e_hists:
        all_hists.append(hist)
    for hist in p_hists:
        all_hists.append(hist)
    for hist in g_hists:
        all_hists.append(hist)

    # Other histograms from full df
    h_helicity = df.Histo1D(("helicity", "Helicity Distribution; Helicity; Counts", 10, -5, 5), "event.helicity")

    h_Q2 = df.Histo1D(("Q2", "Q^{2}; Q^{2} [GeV^2]; Counts", 100, -1, 8), "dis.Q2")
    h_W  = df.Histo1D(("W", "W; W [GeV]; Counts", 200, 0, 10), "dis.W")
    h_Xb = df.Histo1D(("Xb", "x_{B}; x_{B}; Counts", 200, 0, 1), "dis.Xb")
    h_Q2_vs_Xb = df.Histo2D(("Q2_vs_Xb", "Q^{2} vs. x_{B}; x_{B}; Q^{2} [GeV^{2}]", 500, 0, 1, 500, 0, 7), "dis.Xb", "dis.Q2")

    # Open output file for writing
    out_filename = os.path.basename(args.input_file)
    outfile = ROOT.TFile.Open(f"recPlots_{out_filename}", "RECREATE")

    # Write 1D and 2D histograms
    for h in [h_helicity, h_Q2, h_W, h_Xb, h_Q2_vs_Xb]:
        h.Write()
    for hist in all_hists:
        hist.Write()


    outfile.Close()
    print(f"Histograms and canvases saved to rec_{out_filename}")

if __name__ == "__main__":
    main()
