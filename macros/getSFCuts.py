#type: ignore
import ROOT
import argparse
import numpy as np
import matplotlib.pyplot as plt
import json

# Load compiled dictionary for branch variables
ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

parser = argparse.ArgumentParser(description="Compute SF cuts.")
parser.add_argument("input_file", type=str, help="Path to the input ROOT file")
parser.add_argument("--max_events", type=int, help="Maximum number of events to visualize")
parser.add_argument("-q", "--quiet", action="store_true", help="Suppress binned histograms.")
args = parser.parse_args()

file = ROOT.TFile(args.input_file, 'READ')
tree = file.Get("Events")
n_entries = tree.GetEntries()

# Extract electron information
p_list = []
SF_list = []
sector_list = []
pid_list = []

for i in range(n_entries):
    if args.max_events and i >= args.max_events:
        break
    tree.GetEntry(i)
    pid = tree.rec.pid
    if pid != 11:
        continue
    sector = tree.rec.sector
    if sector < 0:
        continue
    p = tree.rec.p
    E_tot = tree.rec.E_PCAL + tree.rec.E_ECIN + tree.rec.E_ECOUT
    if p <= 0 or E_tot <= 0:
        continue
    p_list.append(p)
    SF_list.append(E_tot / p)
    sector_list.append(sector)
    pid_list.append(pid)

# Convert to arrays
p_arr = np.array(p_list)
SF_arr = np.array(SF_list)
sector_arr = np.array(sector_list)
pid_arr = np.array(pid_list)

# Define momentum bins
p_edges = np.linspace(1, 7.5, 23)
p_bin_centers = 0.5 * (p_edges[:-1] + p_edges[1:])

# Prepare storage
mu = [ [] for _ in range(6) ]
sigma = [ [] for _ in range(6) ]

# Loop over sectors and momentum bins
for sec in range(1, 7):
    mask_sector = sector_arr == sec
    SF_sec = SF_arr[mask_sector]
    p_sec = p_arr[mask_sector]

    mu_sec = []
    sigma_sec = []

    for i in range(len(p_edges)-1):
        mask_bin = (p_sec >= p_edges[i]) & (p_sec < p_edges[i+1])
        SF_bin = SF_sec[mask_bin]
        if len(SF_bin) < 10:
            mu_sec.append(np.nan)
            sigma_sec.append(np.nan)
            continue

        # Fit Gaussian
        h_name = f"h_sec{sec}_bin{i}"
        h = ROOT.TH1D(h_name, h_name, 75, SF_bin.min(), SF_bin.max())
        for val in SF_bin:
            h.Fill(val)
        fit = ROOT.TF1("fit", "gaus", SF_bin.min(), SF_bin.max())
        h.Fit(fit, "Q")
        mu_sec.append(fit.GetParameter(1))
        sigma_sec.append(fit.GetParameter(2))

    mu[sec-1] = np.array(mu_sec)
    sigma[sec-1] = np.array(sigma_sec)

# --- Plot 2x3 panel layout with 2D histogram + mu/3sigma bands ---
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
axes = axes.flatten()

fit_results = {}  # store coefficients here

for sec in range(1, 7):
    ax = axes[sec-1]
    mask_sector = (sector_arr == sec)
    SF_sec = SF_arr[mask_sector]
    p_sec = p_arr[mask_sector]

    # 2D histogram as background
    hb = ax.hist2d(p_sec, SF_sec, bins=[75, 75], cmap="cividis", cmin=1)

    # Scatter the raw mu and ±3σ points
    ax.scatter(p_bin_centers, mu[sec-1], color="#E31A1C", s=30, label="mu(p)")
    ax.scatter(p_bin_centers, mu[sec-1]+3*sigma[sec-1], color="#66FFFF", s=20, label="+3σ")
    ax.scatter(p_bin_centers, mu[sec-1]-3*sigma[sec-1], color="#66FFFF", s=20, label="-3σ")

    # Fit mu with a 2nd-degree polynomial (skip NaNs)
    valid_mask = ~np.isnan(mu[sec-1])
    coeff_mu = np.polyfit(p_bin_centers[valid_mask], mu[sec-1][valid_mask], 2)
    poly_mu = np.poly1d(coeff_mu)

    valid_mask_sigma = ~np.isnan(sigma[sec-1])
    coeff_sigma = np.polyfit(p_bin_centers[valid_mask_sigma], sigma[sec-1][valid_mask_sigma], 2)
    poly_sigma = np.poly1d(coeff_sigma)

    # Draw curves
    x_fit = np.linspace(p_bin_centers[0], p_bin_centers[-1], 200)
    ax.plot(x_fit, poly_mu(x_fit), color="#E31A1C", lw=2, linestyle="-")
    ax.plot(x_fit, poly_mu(x_fit)+3*poly_sigma(x_fit), color="#66FFFF", lw=2, linestyle="--")
    ax.plot(x_fit, poly_mu(x_fit)-3*poly_sigma(x_fit), color="#66FFFF", lw=2, linestyle="--")

    ax.set_title(f"Sector {sec}")
    ax.set_xlabel("p [GeV]")
    ax.set_ylabel("Sampling Fraction")
    ax.grid(True)

    # Save coeffs for this sector
    fit_results[f"sector_{sec}"] = {
        "mu_coeffs": coeff_mu.tolist(),
        "sigma_coeffs": coeff_sigma.tolist()
    }

# Write coefficients to JSON
with open("SF_sigma_cut_params.json", "w") as f:
    json.dump(fit_results, f, indent=2)

# Add a single legend and adjust layout
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', fontsize=12)
plt.tight_layout()
plt.subplots_adjust(top=0.92)

# Save figure
plt.savefig("SF_3sigma_per_sector.png", dpi=300)
plt.close()
