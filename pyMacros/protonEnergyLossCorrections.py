#type: ignore
import ROOT
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.offsetbox import AnchoredText
import numpy as np

# Load your compiled dictionary for branch variables
ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

def get_adaptive_edges_from_df(df, varname, n_bins, min_val, max_val, fine_bins=500):
    # Fill fine histogram
    hist = df.Histo1D(("tmp", "tmp", fine_bins, min_val, max_val), varname)
    hist_val = hist.GetValue()
    total_entries = hist_val.GetEntries()
    step = total_entries / n_bins
    
    edges = [hist_val.GetXaxis().GetXmin()]
    cumulative = 0
    target = step
    
    for b in range(1, hist_val.GetNbinsX() + 1):
        cumulative += hist_val.GetBinContent(b)
        if cumulative >= target:
            edge = hist_val.GetBinLowEdge(b + 1)
            if edge > edges[-1]:
                edges.append(edge)
                target += step
    
    if edges[-1] < hist_val.GetXaxis().GetXmax():
        edges.append(hist_val.GetXaxis().GetXmax())
    
    hist_val.Delete()
    return edges

# ---------------------------
# Parse command-line arguments
# ---------------------------
parser = argparse.ArgumentParser(description="Compute proton energy loss corrections from input MC ROOT file.")
parser.add_argument("input_file", type=str, help="Path to the input ROOT file")
parser.add_argument("--max_events", type=int, help="Maximum number of events to visualize")
parser.add_argument("--isCD", action="store_true", help="Proton detected in CD")
parser.add_argument("--doInt", action="store_true", help="Skip binned case, do only the integrated case.")
args = parser.parse_args()
# ---------------------------
# Open ROOT file and create RDataFrame
# ---------------------------
file = ROOT.TFile(args.input_file, 'READ')
tree = file.Get("Events")
df = ROOT.RDataFrame(tree)

# Select protons
radtoDeg = 180 / np.pi
df_pro = df.Filter("rec.pid == 2212 && rec.charge == 1")
df_pro = df_pro.Define("delta_p", "gen.p - rec.p")
df_pro = df_pro.Define("delta_theta", "gen.theta - rec.theta")
df_pro = df_pro.Define("delta_phi", "gen.phi - rec.phi")
df_pro = df_pro.Define("theta_deg", f"rec.theta * {radtoDeg}")
df_pro = df_pro.Define("delta_theta_deg", f"delta_theta * {radtoDeg}")
df_pro = df_pro.Define("delta_phi_deg", f"delta_phi * {radtoDeg}")

if args.max_events is not None:
        df_pro = df_pro.Range(args.max_events)

# ---------------------------
# Integrated case
# ---------------------------
# integrated_hists = []
# canvas = ROOT.TCanvas("c_delta_p_vs_p", "Proton Delta p vs p", 1200, 600)
# canvas.Divide(2, 1)

# df_FD = df_pro.Filter("rec.det == 1")
# hname = "h_delta_p_vs_p_FD"
# hist_FD = df_FD.Histo2D((hname, f"FD Proton #Delta p vs p; p_{{rec}} [GeV]; #Delta p [GeV]",
#                         75, 0, 3, 75, -0.04, 0.2), "rec.p", "delta_p")

# integrated_hists.append(hist_FD)

# canvas.cd(1)
# pad = ROOT.gPad
# pad.SetLeftMargin(0.13)
# pad.SetRightMargin(0.13)
# pad.SetLogz()
# hist_FD.Draw("COLZ")
# hist_FD.SetStats(False) 

# df_CD = df_pro.Filter("rec.det == 2")
# hname = "h_delta_p_vs_p_CD"
# hist_CD = df_CD.Histo2D((hname, f"CD Proton #Delta p vs p; p_{{rec}} [GeV]; #Delta p [GeV]",
#                         75, 0, 3, 75, -0.1, 0.1), "rec.p", "delta_p")

# integrated_hists.append(hist_CD)

# canvas.cd(2)
# pad = ROOT.gPad
# pad.SetLeftMargin(0.13)
# pad.SetRightMargin(0.13)
# pad.SetLogz()
# hist_CD.Draw("COLZ")
# hist_CD.SetStats(False) 

# canvas.SaveAs("proton_delta_p_vs_p_int.png")

# if args.doInt:
#     exit()
# ---------------------------
# Theta binning
# ---------------------------
theta_edges = np.linspace(35, 70, 13) if args.isCD else np.linspace(12, 44, 25)

# if args.isCD:
#     theta_edges = get_adaptive_edges_from_df(df_pro, "theta_deg", n_bins=10, min_val=35, max_val=70)
# else:
#     theta_edges = get_adaptive_edges_from_df(df_pro, "theta_deg", n_bins=24, min_val=6, max_val=44)

# n_bins = len(theta_edges) - 1

def make_theta_binned_fit(
    df,
    variable_x: str,
    variable_y: str,
    theta_edges: list,
    is_cd: bool,
):
    """
    Loop over theta bins and produce 2D histograms with profiles and fits.

    Parameters:
    -----------
    df : ROOT.RDataFrame
        Input dataframe
    variable_x : str
        Name of x variable (e.g., "rec.p")
    variable_y : str
        Name of y variable (e.g., "delta_p")
    theta_edges : list of float
        Bin edges for theta
    is_cd : bool
        Flag for Central Detector (CD) or Forward Detector (FD)

    Returns:
    --------
    theta_centers : list of float
    fit_params : dict of lists
        Contains 'A', 'B', and optionally 'C' values and errors
    objects : list
        List of ROOT histograms, profiles, and fit functions
    """
    n_bins = len(theta_edges) - 1
    n_rows = 2 if is_cd else 4
    n_cols = n_bins // n_rows
    canvas = ROOT.TCanvas("c_thetaBins", "Proton Delta phi vs p in theta bins", 3000, 400 * n_rows)
    canvas.Divide(n_cols, n_rows)
    objects = []
    theta_centers = []
    fit_params = {"A": [], "A_err": [], "B": [], "B_err": []}
    if is_cd:
        fit_params["C"] = []
        fit_params["C_err"] = []

    for i in range(n_bins):
        tmin, tmax = theta_edges[i], theta_edges[i+1]

        # Detector-specific cuts and fit formula
        if is_cd:
            df_cut = df.Filter(f"rec.det == 2 && theta_deg > {tmin} && theta_deg < {tmax}")
            fit_formula = "[0] + [1]*x + [2]*(x*x)"
            y_range = (-0.2, 0.2)
        else:
            df_cut = df.Filter(f"rec.det == 1 && theta_deg > {tmin} && theta_deg < {tmax}")
            fit_formula = "[0] + [1]/x"
            y_range = (-1, 1)

        # 2D histogram
        hname = f"h_{variable_y}_vs_{variable_x}_theta_{int(tmin)}_{int(tmax)}"
        hist = df_cut.Histo2D((
            hname,
            f"{tmin:.1f} < #theta < {tmax:.1f}; p_{{rec}}; #Delta #phi",
            75, 0, 5, 75, y_range[0], y_range[1]
        ), variable_x, variable_y)
        objects.append(hist)

        # Draw on canvas
        canvas.cd(i+1)
        pad = ROOT.gPad
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.13)
        pad.SetLogz()
        hist.Draw("COLZ")

        # Profile and fit
        prof = hist.GetValue().ProfileX()
        prof.SetLineColor(ROOT.kBlack)
        prof.SetMarkerColor(ROOT.kBlack)
        prof.SetMarkerStyle(10)
        objects.append(prof)

        # Fit range
        nbins = prof.GetNbinsX()
        xmin = xmax = None
        for b in range(1, nbins+1):
            if prof.GetBinEntries(b) > 0:
                x = prof.GetBinCenter(b)
                if xmin is None: xmin = x
                xmax = x

        if xmin is not None and xmax is not None and xmax > xmin:
            fit_func = ROOT.TF1(f"f_fit_{i}", fit_formula, xmin, xmax)
            prof.Fit(fit_func, "R")

            # Store fit parameters
            theta_centers.append(0.5*(tmin + tmax))
            fit_params["A"].append(fit_func.GetParameter(0))
            fit_params["A_err"].append(fit_func.GetParError(0))
            fit_params["B"].append(fit_func.GetParameter(1))
            fit_params["B_err"].append(fit_func.GetParError(1))
            if is_cd:
                fit_params["C"].append(fit_func.GetParameter(2))
                fit_params["C_err"].append(fit_func.GetParError(2))

            objects.append(fit_func)

        # Draw overlay
        hist.Draw("COLZ")
        hist.SetStats(False)
        #prof.Draw("SAME")
        #fit_func.Draw("SAME")

    canvas.SaveAs(f"{variable_y}_vs_{variable_x}_thetaBins.png")
    return theta_centers, fit_params, objects

if args.isCD:
    theta_centers, fit_params, objects = make_theta_binned_fit(df_pro, "rec.p", "delta_phi_deg", theta_edges, is_cd=True)
theta_centers, fit_params, objects = make_theta_binned_fit(df_pro, "rec.p", "delta_phi_deg", theta_edges, is_cd=False)


fitA_vals = fit_params["A"]
fitA_errs = fit_params["A_err"]

fitB_vals = fit_params["B"]
fitB_errs = fit_params["B_err"]

# if args.isCD:
#     fitC_vals = fit_params["C"]
#     fitC_errs = fit_params["C_err"]

# ---------------------------
# Plot fit parameters vs theta
# ---------------------------

# Function to setup grid and minor ticks
def setup_axes(ax, x_minor=1.0, y_minor=0.1):
    ax.minorticks_on()
    ax.grid(which='major', linestyle='-', linewidth=0.8)
    ax.grid(which='minor', linestyle='--', linewidth=0.5)
    ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    ax.yaxis.set_minor_locator(MultipleLocator(y_minor))


def safe_polyfit(x, y, y_err=None, order=2):
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    
    if y_err is not None:
        y_err = np.array(y_err, dtype=float)
        mask_valid = (
            np.isfinite(x) & np.isfinite(y) &
            np.isfinite(y_err) & (y_err > 0)
        )
        if mask_valid.sum() < order + 1:
            return None  # not enough valid points to fit
        return np.polyfit(x[mask_valid], y[mask_valid], order, w=1/y_err[mask_valid])
    else:
        mask_valid = np.isfinite(x) & np.isfinite(y)
        if mask_valid.sum() < order + 1:
            return None
        return np.polyfit(x[mask_valid], y[mask_valid], order)


def fit_and_overlay(ax, x, y, y_err=None, order=2, color='red'):
    """
    Fit y(x) to a polynomial of given order, overlay fit on ax, and show fit parameters + reduced chi2.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to plot on
    x : array-like
        x values
    y : array-like
        y values
    y_err : array-like, optional
        Uncertainties on y values; if None, all points treated equally
    order : int
        Polynomial order (default 2)
    color : str
        Color of fit line
    """
    # Fit polynomial with optional weights
    if y_err is not None:
        #weights = 1 / np.array(y_err)
        #p = np.polyfit(x, y, order, w=weights)
        p = safe_polyfit(x, y, y_err, order)
    else:
        #p = np.polyfit(x, y, order)
        p = safe_polyfit(x, y, order=order)
    
    # Fitted y values
    y_fit = np.polyval(p, x)
    
    # Compute reduced chi2
    if y_err is not None:
        chi2 = np.sum(((y - y_fit)/y_err)**2)
    else:
        chi2 = np.sum((y - y_fit)**2)
    ndf = len(x) - len(p)
    reduced_chi2 = chi2 / ndf if ndf > 0 else np.nan
    
    # Overlay fit line
    x_dense = np.linspace(min(x), max(x), 200)
    y_dense = np.polyval(p, x_dense)
    ax.plot(x_dense, y_dense, color=color, linestyle='-', label=f'{order}-order poly fit')
    
    # Stat box
    param_text = '\n'.join([f'p{i} = {pi:.5f}' for i, pi in enumerate(p)])
    param_text += f'\nχ²/ndf = {reduced_chi2:.2f}'
    at = AnchoredText(param_text, loc='upper right', prop=dict(size=9), frameon=True)
    at.patch.set_alpha(0.3)
    ax.add_artist(at)
    
    return p

theta = np.array(theta_centers, dtype=float)
if args.isCD:
    mask = theta >= 40
else:
    mask = theta <= 40
A = np.array(fitA_vals, dtype=float)
A_errs = np.array(fitA_errs, dtype=float)
B = np.array(fitB_vals, dtype=float)
B_errs = np.array(fitB_errs, dtype=float)
fig, axs = plt.subplots(1, 2, figsize=(9, 6)) 
# if args.isCD:
#     C = np.array(fitC_vals, dtype=float)
#     C_errs = np.array(fitC_errs, dtype=float)
#     fig, axs = plt.subplots(1, 3, figsize=(9, 6)) 

# A vs theta
axs[0].errorbar(theta, A, yerr=A_errs, fmt='.')
fitOrder = 2 #if args.isCD else 1
fit_A = fit_and_overlay(axs[0], theta, A, y_err=A_errs, order=fitOrder, color='red')
axs[0].set_ylabel(r"$A_{\phi}$")
axs[0].set_xlabel("θ [deg]")
axs[0].set_title(r"$A_{\phi}(θ)$")
axs[0].grid(True)
if args.isCD:
    axs[0].set_ylim(-1, 1)
else:
    axs[0].set_ylim(-2, 2)
setup_axes(axs[0])

# B vs theta
axs[1].errorbar(theta, B, yerr=B_errs, fmt='.', color='orange')
fit_B = fit_and_overlay(axs[1], theta, B, y_err=B_errs, order=fitOrder, color='red')
axs[1].set_ylabel(r"$B_{\phi}$")
axs[1].set_xlabel("θ [deg]")
axs[1].set_title(r"$B_{\phi}(θ)$")
axs[1].grid(True)
if args.isCD:
    axs[1].set_ylim(-1, 1)
else:
    axs[1].set_ylim(-2, 2)
setup_axes(axs[1])

# if args.isCD:
#     # C vs theta
#     axs[2].errorbar(theta, C, yerr=C_errs, fmt='.', color='green')
#     fit_C = fit_and_overlay(axs[2], theta[mask], C[mask], y_err=C_errs[mask], order=2, color='red')
#     axs[2].set_ylabel(r"$C_{\theta}$")
#     axs[2].set_xlabel("θ [deg]")
#     axs[2].set_title(r"$C_{\theta}(θ)$")
#     axs[2].grid(True)
#     axs[2].set_ylim(-0.3, 0.3)
#     setup_axes(axs[2])

plt.tight_layout()
plt.savefig("proton_delta_theta_vs_p_fitParams.png", dpi=300)
plt.close()

print("All histograms, fits, and parameter plots saved.")
