#type: ignore
import ROOT
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.offsetbox import AnchoredText
import numpy as np
import json 
import re

# Load compiled dictionary for branch variables
ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

def get_adaptive_edges_from_df(df, varname, n_bins, min_val, max_val, fine_bins=10000000):
    # Fine histogram
    hist = df.Histo1D(("tmp", "tmp", fine_bins, min_val, max_val), varname).GetValue()
    total_entries = hist.GetEntries()
    step = total_entries / n_bins

    edges = [min_val]
    cumulative = 0
    target = step

    for b in range(1, hist.GetNbinsX() + 1):
        cumulative += hist.GetBinContent(b)
        while cumulative >= target and len(edges) < n_bins:
            # place edge at the upper edge of this bin
            edge = hist.GetBinLowEdge(b) + hist.GetBinWidth(b)
            if edge > edges[-1]:
                edges.append(edge)
            target += step

    edges.append(max_val)  # always append max
    hist.Delete()
    return edges

def get_adaptive_edges_numpy(df, varname, n_bins, min_val, max_val):
    # grab the column as a numpy array
    x = df.AsNumpy([varname])[varname]
    # cut to desired range
    x = x[(x >= min_val) & (x <= max_val)]
    # quantile edges
    edges = np.quantile(x, np.linspace(0, 1, n_bins+1))
    return edges.tolist()

def plot_integrated_proton_deltap_vs_p(df_pro):
    """
    Plot integrated proton Δp vs p for FD and CD detectors.

    df_pro : ROOT.RDataFrame
        DataFrame containing proton reconstruction and delta_p variables.
    """
    integrated_hists = []

    # Create canvas
    canvas = ROOT.TCanvas("c_delta_p_vs_p", "Proton Delta p vs p", 1200, 600)
    canvas.Divide(2, 1)

    # ---------------------------
    # Forward Detector (FD)
    # ---------------------------
    df_FD = df_pro.Filter("rec.det == 1")
    hist_FD = df_FD.Histo2D(
        ("h_delta_p_vs_p_FD",
         "FD Proton #Delta p vs p; p_{rec} [GeV]; #Delta p [GeV]",
         75, 0, 3, 75, -0.04, 0.2),
        "rec.p", "delta_p"
    )
    integrated_hists.append(hist_FD)

    canvas.cd(1)
    pad = ROOT.gPad
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.13)
    pad.SetLogz()
    hist_FD.Draw("COLZ")
    hist_FD.SetStats(False)

    # ---------------------------
    # Central Detector (CD)
    # ---------------------------
    df_CD = df_pro.Filter("rec.det == 2")
    hist_CD = df_CD.Histo2D(
        ("h_delta_p_vs_p_CD",
         "CD Proton #Delta p vs p; p_{rec} [GeV]; #Delta p [GeV]",
         75, 0, 3, 75, -0.1, 0.1),
        "rec.p", "delta_p"
    )
    integrated_hists.append(hist_CD)

    canvas.cd(2)
    pad = ROOT.gPad
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.13)
    pad.SetLogz()
    hist_CD.Draw("COLZ")
    hist_CD.SetStats(False)

    # Save the plot
    canvas.SaveAs("proton_delta_p_vs_p_int.png")

    return integrated_hists

def theta_binned_fit(df, theta_edges: list, variable_y: str, fit_formula: str, y_range: tuple, det: int, nrows=4, quiet=False):
    """
    Loop over theta bins, produce 2D histograms, profiles, and fits.

    Parameters
    ----------
    df : ROOT.RDataFrame
        Input dataframe
    variable_x : str
        X variable (e.g., "rec.p")
    variable_y : str
        Y variable (e.g., "delta_p")
    theta_edges : list[float]
        Bin edges for theta
    fit_formula : str
        ROOT TF1 formula string (e.g. "[0] + [1]*x", "[0] + [1]/x + [2]*x*x")
    y_range : tuple[float, float]
        Y-axis range for histograms
    det_id : int
        Detector ID to filter (1=FD, 2=CD)
    canvas_nrows: int
    
    Returns
    -------
    theta_centers : list[float]
        Midpoints of theta bins
    fit_params : dict
        Keys: "p0", "p1", ..., "pN" and corresponding "_err" lists
    objects : list
        ROOT objects created (histograms, profiles, fits)
    """
    n_bins = len(theta_edges) - 1

    # Determine number of fit parameters from formula
    param_indices = sorted(set(int(n) for n in re.findall(r"\[(\d+)\]", fit_formula)))
    n_params = max(param_indices) + 1

    fit_params = {}
    for i in range(n_params):
        fit_params[f"p{i}"] = []
        fit_params[f"p{i}_err"] = []
    
    objects = []
    theta_centers = []

    if not quiet:
        canvas = ROOT.TCanvas("c_thetaBins", "theta binned fits", 3200, 400 * nrows)
        ncols = n_bins // nrows
        canvas.Divide(ncols, nrows)

    for i in range(n_bins):
        tmin, tmax = theta_edges[i], theta_edges[i+1]
        df_cut = df.Filter(f"rec.det == {det} && theta_deg > {tmin} && theta_deg < {tmax}")

        # 2D histogram
        hname = f"h_{variable_y}_vs_p_theta_{int(tmin)}_{int(tmax)}"
        hist = df_cut.Histo2D((
            hname,
            f"{tmin:.1f} < #theta < {tmax:.1f}; p_{{rec}}; #Delta {root_symbol_map.get(variable_y, variable_y)}",
            75, 0, 5, 75, y_range[0], y_range[1]
        ), "rec.p", variable_y)
        hist.GetXaxis().SetTitleSize(0.06)  # X-axis title
        hist.GetYaxis().SetTitleSize(0.06)  # Y-axis title
        objects.append(hist)

        if not quiet:
            canvas.cd(i+1)
            pad = ROOT.gPad
            pad.SetLeftMargin(0.135)
            pad.SetRightMargin(0.135)
            pad.SetBottomMargin(0.13)
            pad.SetLogz()
            hist.Draw("COLZ")

        # Profile and fit
        prof = hist.GetValue().ProfileX()
        prof.SetLineColor(ROOT.kBlack)
        prof.SetMarkerColor(ROOT.kBlack)
        prof.SetMarkerStyle(10)
        objects.append(prof)

        nbins = prof.GetNbinsX()
        xmin, xmax = None, None
        for b in range(1, nbins+1):
            if prof.GetBinEntries(b) > 0:
                x = prof.GetBinCenter(b)
                if xmin is None: xmin = x
                xmax = x

        if xmin is not None and xmax is not None and xmax > xmin:
            fit_func = ROOT.TF1(f"f_fit_{i}", fit_formula, xmin, xmax)
            prof.Fit(fit_func, "RQ")  
            theta_centers.append(0.5 * (tmin + tmax))

            # Store parameters dynamically
            for p in range(n_params):
                fit_params[f"p{p}"].append(fit_func.GetParameter(p))
                fit_params[f"p{p}_err"].append(fit_func.GetParError(p))

            objects.append(fit_func)
        
        if not quiet:
            # Draw overlay
            hist.Draw("COLZ")
            hist.SetStats(False)
            prof.Draw("SAME")
            fit_func.Draw("SAME")

    if not quiet:
        detector = "FD" if det == 1 else "CD"
        canvas.SaveAs(f"{detector}_{variable_y}_vs_p_thetaBins.png")
        
    return theta_centers, fit_params, objects

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

def fit_and_overlay(ax, x, y, y_err=None, order=2, color='red', theta_threshold=None, color_below='blue', color_above='green'):
    """
    Fit y(x) to a polynomial of given order, overlay fit on ax, and show fit parameters + reduced chi2.
    Can optionally split fits below and above a theta_threshold.

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
        Color of fit line if no threshold split
    theta_threshold : float, optional
        Split x values at this threshold and fit separately
    """
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    if y_err is not None:
        y_err = np.array(y_err, dtype=float)
    
    all_coeffs = {}

    # Helper function for fitting and annotating a single region
    def _fit_region(x_region, y_region, y_err_region, region_key, line_color, box_loc='upper right', y_offset=0.0):
        p = safe_polyfit(x_region, y_region, y_err_region, order)
        if p is None:
            return None
        y_fit = np.polyval(p, x_region)
        chi2 = np.sum(((y_region - y_fit)/y_err_region)**2) if y_err_region is not None else np.sum((y_region - y_fit)**2)
        ndf = len(x_region) - len(p)
        reduced_chi2 = chi2 / ndf if ndf > 0 else np.nan
        x_dense = np.linspace(min(x_region), max(x_region), 200)
        y_dense = np.polyval(p, x_dense)
        ax.plot(x_dense, y_dense, color=line_color, linestyle='-', label=f'Fit {line_color}')
        # Stat box
        param_text = '\n'.join([f'p{i} = {pi:.5f}' for i, pi in enumerate(p[::-1])])
        param_text += f'\nχ²/ndf = {reduced_chi2:.2f}'
        at = AnchoredText(param_text, loc=box_loc, prop=dict(size=9), frameon=True)
        at.patch.set_alpha(0.3)

        # Shift the box vertically by a fraction of the axis height
        at.set_bbox_to_anchor((1, 1 - y_offset), transform=ax.transAxes)
        ax.add_artist(at)
        all_coeffs[region_key] = p[::-1].tolist()
        return p[::-1]

    if theta_threshold is None:
        _fit_region(x, y, y_err, "all", color)
    else:
        mask_below = x <= theta_threshold
        mask_above = x > theta_threshold
        if np.any(mask_below):
            _fit_region(x[mask_below], y[mask_below], y_err[mask_below] if y_err is not None else None, region_key="below", line_color=color_below, y_offset=0)
        if np.any(mask_above):
            _fit_region(x[mask_above], y[mask_above], y_err[mask_above] if y_err is not None else None, region_key="above", line_color=color_above, y_offset=0.25)

    return all_coeffs

# ---------------------------
# MAIN WORKFLOW:
# ---------------------------
parser = argparse.ArgumentParser(description="Compute proton energy loss corrections from input MC ROOT file.")
parser.add_argument("input_file", type=str, help="Path to the input ROOT file")
parser.add_argument("--max_events", type=int, help="Maximum number of events to visualize")
parser.add_argument("-q", "--quiet", action="store_true", help="Suppress binned histograms.")
parser.add_argument("--isCD", action="store_true", help="Proton detected in CD")
args = parser.parse_args()
 
file = ROOT.TFile(args.input_file, 'READ')
tree = file.Get("Events")
df = ROOT.RDataFrame(tree)

radtoDeg = 180 / np.pi
# Select protons:
df_pro = df.Filter("rec.pid == 2212 && rec.charge == 1")
df_pro = df_pro.Define("theta_deg", f"rec.theta * {radtoDeg}")
df_pro = df_pro.Define("delta_p", "gen.p - rec.p")
df_pro = df_pro.Define("delta_theta", f"(gen.theta - rec.theta) * {radtoDeg}")
df_pro = df_pro.Define("delta_phi", f"TMath::ATan2(TMath::Sin(gen.phi - rec.phi), TMath::Cos(gen.phi - rec.phi)) * {radtoDeg}")

mpl_symbol_map  = {"delta_p": "p", "delta_theta": r"\theta", "delta_phi": r"\phi"}
root_symbol_map = {"delta_p": "p", "delta_theta": "#theta", "delta_phi": "#phi"}

if args.max_events is not None:
    df_pro = df_pro.Range(args.max_events)

plot_integrated_proton_deltap_vs_p(df_pro)

# ---------------------------
# Theta binning
# ---------------------------
#theta_edges = np.linspace(35, 70, 13) if args.isCD else np.linspace(12, 44, 25)

if args.isCD:
    theta_edges = get_adaptive_edges_numpy(df_pro, "theta_deg", n_bins=10, min_val=40, max_val=58)
else:
    theta_edges = get_adaptive_edges_numpy(df_pro, "theta_deg", n_bins=24, min_val=15, max_val=40)

# ---------------------------
# Configs for binned fits
# ---------------------------
fit_configs = {
    "delta_p": ["[0] + [1]/x + [2]/(x*x)", (-0.06, 0.06), 1, 4],
    "delta_theta": ["[0] + [1]/x", (-1, 1), 1, 4],
    "delta_phi": ["[0] + [1]/x + [2]/(x*x)", (-1, 1), 1, 4],
}

if args.isCD:
    fit_configs["delta_p"] = ["[0] + [1]*x + [2]*(x*x)", (-0.2, 0.2), 2, 2]
    fit_configs["delta_theta"] = ["[0] + [1]/x", (-0.2, 0.2), 2, 2]
    fit_configs["delta_phi"] = ["[0] + [1]/x + [2]/(x*x)", (-0.2, 0.2), 2, 2]
    
# ---------------------------
# Run all theta-binned fits
# ---------------------------
results = {}  # store bin centers & fit info
for var, args_list in fit_configs.items():
    centers, params, objs = theta_binned_fit(df_pro, theta_edges, var, *args_list, args.quiet)
    results[var] = {
        "centers": np.array(centers, dtype=float),
        "params": params,
        "objects": objs,
    }

# ---------------------------
# Plot fit parameters vs theta
# ---------------------------

# Define polynomial fit orders for each (variable, detector) combination
fit_orders = {
    # variable : { detector : order }
    "delta_p": {"FD": 1, "CD": 2},
    "delta_theta": {"FD": 3, "CD": 2},
    "delta_phi": {"FD": 4, "CD": 2},
}
detector_key = "CD" if args.isCD else "FD"

THRESHOLD_THETA = None #deg

alphanum_map = {"p0": "A", "p1": "B", "p2": "C"}

fit_results_export = {}

for var, res in results.items():
    centers = res["centers"]
    params = res["params"]

    # Auto-detect number of fit params
    param_names = [k for k in params.keys() if not k.endswith("_err")]
    fig, axs = plt.subplots(1, len(param_names), figsize=(5*len(param_names), 5))

    if len(param_names) == 1:
        axs = [axs]

    for ax, pname in zip(axs, param_names):
        y = np.array(params[pname], dtype=float)
        y_err = np.array(params[f"{pname}_err"], dtype=float)

        ax.errorbar(centers, y, yerr=y_err, fmt='.', label='Theta-binned params')
         # Lookup fit order — fallback to 1 if not defined
        fit_order = fit_orders.get(var, {}).get(detector_key, 1)
        coeffs = fit_and_overlay(ax, centers, y, y_err, order=fit_order, color='red', theta_threshold=THRESHOLD_THETA)

        display_name = alphanum_map.get(pname, pname)
        key = fr"${display_name}_{{{mpl_symbol_map.get(var)}}}$"   
        fit_results_export[key] = coeffs
        ax.set_ylabel(rf"${display_name}_{{{mpl_symbol_map.get(var)}}}$")
        ax.set_xlabel(r"$\theta$ [deg]")
        ax.set_title(fr"${display_name}_{{{mpl_symbol_map.get(var)}}}(\theta)$")
        setup_axes(ax)

    plt.tight_layout()
    plt.savefig(f"{detector_key}_{var}_fitParams_vs_theta.png", dpi=300)
    plt.close()

# Write to JSON file
with open(f"{detector_key}_proton_kinCorr_coeffs.json", "w") as f:
    json.dump(fit_results_export, f, indent=2)

print("Saved correction coefficients to proton_kinCorr_coeffs.json")
print("All histograms, fits, and parameter plots saved.")
