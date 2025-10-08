#!/usr/bin/env python3
# type:ignore
import os
import ROOT
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

if len(sys.argv) < 4:
    print("Usage: python quilt2D.py <input.root> <t_index> <chi2_threshold>")
    sys.exit(1)

input_file = sys.argv[1]
chi2_threshold = float(sys.argv[2])
t_index = int(sys.argv[3])

root_file = ROOT.TFile.Open(input_file)
if not root_file or root_file.IsZombie():
    print(f"Could not open {input_file}")
    sys.exit(1)

# --- detect max q and xb bins by scanning ---
q_bins, xb_bins = set(), set()
for key in root_file.GetListOfKeys():
    name = key.GetName()
    if name.startswith("gr_phi_q") and f"_t{t_index}" in name:
        try:
            q_str = name.split("_")[2]  # q0
            xb_str = name.split("_")[3] # xb0
            q_bins.add(int(q_str.replace("q","")))
            xb_bins.add(int(xb_str.replace("xb","")))
        except Exception:
            pass

if not q_bins or not xb_bins:
    print(f"No graphs found for t index {t_index}")
    sys.exit(0)

q_bins = sorted(list(q_bins))
xb_bins = sorted(list(xb_bins))

print(f"Found Q bins: {q_bins}")
print(f"Found xB bins: {xb_bins}")

nq = len(q_bins)
nxb = len(xb_bins)

fig, axes = plt.subplots(nq, nxb, figsize=(5*nxb, 3*nq), squeeze=False)

for i, q in enumerate(q_bins):
    for j, xb in enumerate(xb_bins):
        # invert Q² axis to have increasing Q² upward
        row = nq - 1 - i
        col = j
        ax = axes[row, col]
        graph_name = f"gr_phi_q{q}_xb{xb}_t{t_index}"
        gr = root_file.Get(graph_name)
        if not gr:
            ax.set_visible(False)
            print(f"Missing: {graph_name}")
            continue

        # --- extract points using GetX()/GetY() ---
        n = gr.GetN()
        x = np.array([gr.GetX()[k] for k in range(n)])
        y = np.array([gr.GetY()[k] for k in range(n)])
        yerr = np.array([gr.GetErrorY(k) for k in range(n)])

        # --- scale down by 1e6 ---
        y *= 1e-6
        yerr *= 1e-6

        # Fit A + B cos φ + C cos 2φ
        chi2_ndf = np.nan
        try:
            def model(phi, A, B, C):
                return A + B*np.cos(np.radians(phi)) + C*np.cos(2*np.radians(phi))
            mask = yerr > 0
            if np.any(mask):
                pars, cov = curve_fit(model, x[mask], y[mask], sigma=yerr[mask], absolute_sigma=True,  bounds=([0, -np.inf, -np.inf], [np.inf, np.inf, np.inf]))
                chi2 = np.sum(((y[mask] - model(x[mask], *pars))/yerr[mask])**2)
                ndf = np.sum(mask) - len(pars)
                chi2_ndf = chi2/ndf if ndf>0 else np.nan
            else:
                pars, cov = None, None
        except Exception:
            pars, cov = None, None

        # Skip panels with too-high chi2/ndf
        if chi2_ndf is not np.nan and chi2_ndf > chi2_threshold:
            ax.set_visible(False)
            print(f"Skipping {graph_name} due to chi2/ndf = {chi2_ndf:.2f}")
            continue

        # --- Plot data & fit ---
        ax.errorbar(x, y, yerr=yerr, fmt="o", ms=5, color="black")
        if pars is not None and cov is not None and not (chi2_ndf > chi2_threshold):
            phi_fine = np.linspace(0, 360, 500)
            y_fit = model(phi_fine, *pars)
            y_fit_err = np.sqrt(np.diag(cov)[0]) if cov is not None else 0
            ax.plot(phi_fine, y_fit, "m-")
            ax.fill_between(phi_fine, y_fit-y_fit_err, y_fit+y_fit_err, color="m", alpha=0.2)
            # ax.text(0.95, 0.9, f"$\\chi^2$/ndf={chi2_ndf:.2f}", transform=ax.transAxes,
            #         ha="right", va="top", fontsize=8)

        ax.set_xlim(0, 360)
        ax.set_title(f"Q{q}, xB{xb}", fontsize=9)
        if row == nq-1:
            ax.set_xlabel(r"$\phi$ [deg]")
        if col == 0:
            ax.set_ylabel(r"$d^2\sigma/dtd\phi$")

plt.tight_layout()

# --- Save ---
output_dir = "quilt2D_plots"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
outname = os.path.join(output_dir, f"phi_quilt_t{t_index}.png")
plt.savefig(outname, dpi=300)
print(f"Saved quilt to {outname}")

plt.show()
root_file.Close()



plt.show()
root_file.Close()
