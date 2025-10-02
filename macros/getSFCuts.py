# type: ignore
import ROOT
import argparse
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit

# Load compiled dictionary for branch variables
ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

parser = argparse.ArgumentParser(description="Compute SF cuts.")
parser.add_argument("input_file", type=str, help="Path to the input ROOT file")
parser.add_argument("--gemc", action="store_true", help="Treat input as MC (sector-independent)")
parser.add_argument("--max_events", type=int, help="Maximum number of events to visualize")
parser.add_argument("-q", "--quiet", action="store_true", help="Suppress plotting.")
args = parser.parse_args()

# --- Load TTree --- #
file = ROOT.TFile(args.input_file, 'READ')
tree = file.Get("Events")
n_entries = tree.GetEntries()

# --- Extract electron information --- #
p_list, SF_list, sector_list = [], [], []

for i in range(n_entries):
    if args.max_events and i >= args.max_events:
        break
    tree.GetEntry(i)
    pid = tree.rec.pid
    if pid != 11:  # select electrons
        continue
    p = tree.rec.p
    E_tot = tree.rec.E_PCAL + tree.rec.E_ECIN + tree.rec.E_ECOUT
    if p <= 0 or E_tot <= 0:
        continue
    p_list.append(p)
    SF_list.append(E_tot / p)

    if not args.gemc:
        sector_val = tree.rec.sector
        if sector_val >= 0:
            sector_list.append(sector_val)

# --- Convert to arrays --- #
p_arr = np.array(p_list)
SF_arr = np.array(SF_list)

if not args.gemc:
    sector_arr = np.array(sector_list)

# --- Define momentum bins --- #
p_edges = np.linspace(1, 10, 23)
p_bin_centers = 0.5 * (p_edges[:-1] + p_edges[1:])

for i in range(len(p_edges)-1):
    mask = (p_arr >= p_edges[i]) & (p_arr < p_edges[i+1])
    print(f"Bin {i}, entries = {mask.sum()}")

# --- Fit helper for mu and sigma --- #
def fit_mu_sigma(SF_vals, p_vals):
    mu = []
    sigma = []
    for i in range(len(p_edges)-1):
        mask = (p_vals >= p_edges[i]) & (p_vals < p_edges[i+1])
        SF_bin = SF_vals[mask]
        if len(SF_bin) < 10:
            mu.append(np.nan)
            sigma.append(np.nan)
            continue
        h = ROOT.TH1D("h", "h", 75, SF_bin.min(), SF_bin.max())
        h.SetDirectory(0)
        for val in SF_bin:
            h.Fill(val)
        f = ROOT.TF1("f", "gaus", SF_bin.min(), SF_bin.max())
        h.Fit(f, "Q")
        mu.append(f.GetParameter(1))
        sigma.append(f.GetParameter(2))
    return np.array(mu), np.array(sigma)

# --- Quadratic fit helper with curve_fit --- #
def quadratic(x, a, b, c):
    return a*x**2 + b*x + c

def fit_quadratic(x, y, label=""):
    # Remove NaNs
    mask = ~np.isnan(y)
    x_fit = x[mask]
    y_fit = y[mask]
    if len(x_fit) < 3:
        print(f"[WARN] Not enough points to fit quadratic for {label}")
        return [np.nan, np.nan, np.nan]

    # Initial guesses: a=0, b=0, c=mean(y)
    p0 = [0.0, 0.0, np.nanmean(y_fit)]
    try:
        popt, pcov = curve_fit(quadratic, x_fit, y_fit, p0=p0)
        #print(f"[DEBUG] {label} fit succeeded: coeffs = {popt}")
        return popt.tolist()
    except Exception as e:
        print(f"[ERROR] {label} fit failed: {e}")
        return [np.nan, np.nan, np.nan]

# --- Fit mu and sigma --- #
fit_results = {}
if args.gemc:
    mu, sigma = fit_mu_sigma(SF_arr, p_arr)
    mu_per_sector = [mu for _ in range(6)]
    sigma_per_sector = [sigma for _ in range(6)]
    for sec in range(1, 7):
        coeff_mu = fit_quadratic(p_bin_centers, mu, f"mu sector {sec}")
        coeff_sigma = fit_quadratic(p_bin_centers, sigma, f"sigma sector {sec}")
        fit_results[f"sector_{sec}"] = {"mu_coeffs": coeff_mu, "sigma_coeffs": coeff_sigma}
else:
    mu_per_sector, sigma_per_sector = [], []
    for sec in range(1, 7):
        mask_sec = sector_arr == sec
        mu_sec, sigma_sec = fit_mu_sigma(SF_arr[mask_sec], p_arr[mask_sec])
        mu_per_sector.append(mu_sec)
        sigma_per_sector.append(sigma_sec)

        coeff_mu = fit_quadratic(p_bin_centers, mu_sec, f"mu sector {sec}")
        coeff_sigma = fit_quadratic(p_bin_centers, sigma_sec, f"sigma sector {sec}")
        fit_results[f"sector_{sec}"] = {"mu_coeffs": coeff_mu, "sigma_coeffs": coeff_sigma}

# --- Save JSON --- #
outfile = "SF_sigma_cut_params_GEMC.json" if args.gemc else "SF_sigma_cut_params_REC.json"
with open(outfile, "w") as f:
    json.dump(fit_results, f, indent=2)
print(f"Saved SF coefficients to {outfile}")

# --- Optional plotting --- #
if not args.quiet:
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
    axes = axes.flatten()

    for sec in range(6):
        ax = axes[sec]
        if args.gemc:
            mask_sector = np.ones_like(SF_arr, dtype=bool)
        else:
            mask_sector = (sector_arr == (sec+1))
        SF_sec = SF_arr[mask_sector]
        p_sec = p_arr[mask_sector]

        hb = ax.hist2d(p_sec, SF_sec, bins=[75, 75], cmap="cividis", cmin=1)

        ax.scatter(p_bin_centers, mu_per_sector[sec], color="#E31A1C", s=30, label="mu(p)")
        ax.scatter(p_bin_centers, mu_per_sector[sec]+3*sigma_per_sector[sec], color="#66FFFF", s=20, label="+3σ")
        ax.scatter(p_bin_centers, mu_per_sector[sec]-3*sigma_per_sector[sec], color="#66FFFF", s=20, label="-3σ")

        coeff_mu = fit_results[f"sector_{sec+1}"]["mu_coeffs"]
        coeff_sigma = fit_results[f"sector_{sec+1}"]["sigma_coeffs"]

        if not any(np.isnan(coeff_mu)):
            poly_mu = np.poly1d(coeff_mu)
        else:
            poly_mu = None
        if not any(np.isnan(coeff_sigma)):
            poly_sigma = np.poly1d(coeff_sigma)
        else:
            poly_sigma = None

        x_fit = np.linspace(p_bin_centers[0], p_bin_centers[-1], 200)
        if poly_mu is not None and poly_sigma is not None:
            ax.plot(x_fit, poly_mu(x_fit), color="#E31A1C", lw=2)
            ax.plot(x_fit, poly_mu(x_fit)+3*poly_sigma(x_fit), color="#66FFFF", lw=2, linestyle="--")
            ax.plot(x_fit, poly_mu(x_fit)-3*poly_sigma(x_fit), color="#66FFFF", lw=2, linestyle="--")

        ax.set_title(f"Sector {sec+1}")
        ax.set_xlabel("p [GeV]")
        ax.set_ylabel("Sampling Fraction")
        ax.grid(True)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', fontsize=12)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig("SF_3sigma_per_sector.png", dpi=300)
    plt.close()
