# type: ignore
import argparse
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from datetime import datetime

def current_timestamp():
    now = datetime.now()
    return now.strftime("%m%d_%H%M")

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

###### MAIN BODY STARTS HERE #######

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Compute ALU from input ROOT file.")
parser.add_argument("input_file", type=str, help="Path to the input ROOT file")
parser.add_argument("--isRGA", action="store_true", help="Apply helicity flip for RGA convention")
args = parser.parse_args()

# Open the ROOT file
file = ROOT.TFile(args.input_file, 'READ')

# Access the trees inside the file
tree = file.Get('Events') 

df = ROOT.RDataFrame(tree)

#print(df.GetColumnNames())

# Get adaptive bin edges for each variable
Q2_edges = get_adaptive_edges(tree, "dis.Q2",   n_bins=4, min_val=1.0, max_val=6.5)
Xb_edges = get_adaptive_edges(tree, "dis.Xb",   n_bins=4, min_val=0.1, max_val=0.7)
t_edges  = get_adaptive_edges(tree, "eppi0.t",  n_bins=1, min_val=0.0, max_val=4)

# phi binning: uniform (for now)
phi_edges = np.linspace(0, 360, 13)  # 12 bins, e.g.

print("Adaptive Q2 edges:", Q2_edges)
print("Adaptive Xb edges:", Xb_edges)
print("Adaptive t edges:", t_edges)

nQ2  = len(Q2_edges) - 1
nXb  = len(Xb_edges) - 1
nt   = len(t_edges) - 1
nphi = len(phi_edges) - 1

P = 0.8
if args.isRGA:
    P *= -1

# Prepare arrays for counts: N+ and N- for all bins
N_plus = np.zeros((nQ2, nXb, nt, nphi), dtype=int)
N_minus = np.zeros_like(N_plus)

# Snapshot needed columns into numpy arrays for speed
cols = ["dis.Q2", "dis.Xb", "eppi0.t", "eppi0.trentoPhi", "event.helicity"]
df_numpy = df.AsNumpy(columns=cols)

Q2_vals = df_numpy["dis.Q2"]
Xb_vals = df_numpy["dis.Xb"]
t_vals = df_numpy["eppi0.t"]
phi_vals = (df_numpy["eppi0.trentoPhi"] + 2*np.pi) % (2*np.pi) * 180.0 / np.pi
helicities = df_numpy["event.helicity"]

# Helper function to find bin index
def find_bin(value, edges):
    idx = np.searchsorted(edges, value) - 1
    if idx < 0 or idx >= len(edges) - 1:
        return -1
    return idx

# Fill counts per event
for i in range(len(Q2_vals)):
    iQ2  = find_bin(Q2_vals[i], Q2_edges)
    iXb  = find_bin(Xb_vals[i], Xb_edges)
    it   = find_bin(t_vals[i], t_edges)
    iphi = find_bin(phi_vals[i], phi_edges)

    if -1 in (iQ2, iXb, it, iphi): 
        continue
    h = helicities[i]
    if h == 1:
        N_plus[iQ2, iXb, it, iphi] += 1
    elif h == -1:
        N_minus[iQ2, iXb, it, iphi] += 1

# Now compute asymmetry ALU = (N+ - N-) / (N+ + N-) / P and errors
ALU = np.zeros_like(N_plus, dtype=float)
ALU_err = np.zeros_like(N_plus, dtype=float)

for idx, _ in np.ndenumerate(N_plus):
    Np = N_plus[idx]
    Nm = N_minus[idx]
    Ntot = Np + Nm
    if Ntot > 0:
        ALU[idx] = (Np - Nm) / Ntot / P
        ALU_err[idx] = (2.0 / P) * np.sqrt((Np * Nm) / (Ntot ** 3))
    else:
        ALU[idx] = 0
        ALU_err[idx] = 0

# plot integrated asymmetry over Q2,Xb,t vs phi (sum over first 3 dims)
diff_int = np.sum(N_plus, axis=(0,1,2)) - np.sum(N_minus, axis=(0,1,2))
sum_int  = np.sum(N_plus, axis=(0,1,2)) + np.sum(N_minus, axis=(0,1,2))
with np.errstate(divide='ignore', invalid='ignore'):
    ALU_int = diff_int / sum_int / P
    ALU_int_err = (2.0 / abs(P)) * np.sqrt((np.sum(N_plus, axis=(0,1,2)) * np.sum(N_minus, axis=(0,1,2))) / (sum_int ** 3))

phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])

# Sine function, assume no phase offset: A * sin(phi), phi in degrees
def sin_fit(phi_deg, A):
    phi_rad = np.radians(phi_deg)
    return A * np.sin(phi_rad)

# Mask valid entries (not NaN or inf)
valid = np.isfinite(ALU_int) & np.isfinite(ALU_int_err)
phi_fit = phi_centers[valid]
ALU_fit = ALU_int[valid]
ALU_err = ALU_int_err[valid]

# Perform fit with initial guess: A=0.1
popt, pcov = curve_fit(sin_fit, phi_fit, ALU_fit, sigma=ALU_err, absolute_sigma=True, p0=[0.1])
A_fit = popt[0]
A_fit_err = np.sqrt(np.diag(pcov))[0]

# Plot fit curve
phi_plot = np.linspace(0, 360, 720)
fit_curve = sin_fit(phi_plot, A_fit)

plt.errorbar(phi_centers, ALU_int, yerr=ALU_int_err, fmt='o', color='k', label="RGA fa18t+1 DVPi0P Skim")
plt.plot(phi_plot, fit_curve, linestyle='--', color='r', linewidth=1, label=fr"Fit: $A \sin(\phi)$" + "\n" +
         fr"$A$ = {float(A_fit):.3f} ± {float(A_fit_err):.3f}")
plt.xlabel("φ [deg]")
plt.ylabel(r"$A_{LU}(\phi)$")
plt.title("Integrated Beam Spin Asymmetry")
plt.grid(True)
plt.legend()
plt.tight_layout()
figTitle = f"ALU_phi_fit_{current_timestamp()}.png"
plt.savefig(figTitle, dpi=300)
plt.show()
print(f"integrated ALU vs Phi saved as {figTitle}.")

THRESHOLD = 100  # Minimum total events required to plot

fig, axes = plt.subplots(
    nQ2, nXb,
    figsize=(4 * nXb, 4 * nQ2),
    sharex=True, sharey=True,
    constrained_layout=True
)
fig.suptitle("BSA vs φ in (Q², Xb) bins", fontsize=16)

leftmost_visible_col = [None] * nQ2
bottommost_visible_row = [None] * nXb

visible_axes = []  # track axes that remain

# First pass: determine visibility and record leftmost/bottommost visible indices
for iQ2 in range(nQ2):
    for iXb in range(nXb):
        ax = axes[nQ2 - 1 - iQ2, iXb] if nQ2 > 1 and nXb > 1 else axes[max(iQ2, iXb)]

        Np = np.sum(N_plus[iQ2, iXb, :, :], axis=0)
        Nm = np.sum(N_minus[iQ2, iXb, :, :], axis=0)
        Ntot = Np + Nm

        if np.sum(Ntot) < THRESHOLD:
            fig.delaxes(ax)  # Remove this axis from the figure
            continue

        visible_axes.append((iQ2, iXb, ax))

        if leftmost_visible_col[iQ2] is None or iXb < leftmost_visible_col[iQ2]:
            leftmost_visible_col[iQ2] = iXb
        if bottommost_visible_row[iXb] is None or iQ2 < bottommost_visible_row[iXb]:
            bottommost_visible_row[iXb] = iQ2

# Second pass: plot and set labels, legend, and event counts
for iQ2, iXb, ax in visible_axes:
        ax = axes[nQ2 - 1 - iQ2, iXb] if nQ2 > 1 and nXb > 1 else axes[max(iQ2, iXb)]

        Np = np.sum(N_plus[iQ2, iXb, :, :], axis=0)
        Nm = np.sum(N_minus[iQ2, iXb, :, :], axis=0)
        Ntot = Np + Nm

        with np.errstate(divide='ignore', invalid='ignore'):
            ALU_proj = (Np - Nm) / Ntot / P
            ALU_proj_err = (2.0 / abs(P)) * np.sqrt((Np * Nm) / (Ntot ** 3))
            ALU_proj[~np.isfinite(ALU_proj)] = 0
            ALU_proj_err[~np.isfinite(ALU_proj_err)] = 0

        # Plot data points
        points = ax.errorbar(phi_centers, ALU_proj, yerr=ALU_proj_err, fmt='o', color='k', markersize=4, label='RGA fa18t+1 DVPi0P Skim')

        # Fit curve (optional)
        fit_line = None
        try:
            # Note: Ntot floor set arbitrarily here currently.
            valid = (Ntot > 10)
            popt, _ = curve_fit(
                sin_fit, phi_centers[valid], ALU_proj[valid],
                sigma=ALU_proj_err[valid], absolute_sigma=True, p0=[0.1]
            )
            phi_fit_line = sin_fit(phi_plot, *popt)
            fit_line, = ax.plot(phi_plot, phi_fit_line, linestyle='--', color='r', linewidth=1, label=f'Fit: A={popt[0]:.3f}±{np.sqrt(pcov[0,0]):.3f}')
        except Exception:
            pass

        # Number of events text (top right corner)
        total_events = int(np.sum(Ntot))
        ax.text(
            0.95, 0.90, f"Events: {total_events}",
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes,
            fontsize=9,
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1)
        )

        ax.set_title(
            f"Q² [{Q2_edges[iQ2]:.2f}, {Q2_edges[iQ2+1]:.2f}]\n"
            f"Xb [{Xb_edges[iXb]:.3f}, {Xb_edges[iXb+1]:.3f}]",
            fontsize=10
        )

        if iXb == leftmost_visible_col[iQ2]:
            ax.set_ylabel(r"$A_{LU}(\phi)$")

        if iQ2 == bottommost_visible_row[iXb]:
            ax.set_xlabel("φ [deg]")
            ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
        else:
            ax.tick_params(axis='x', which='both', labelbottom=False)

        ax.set_xlim(0, 360)
        ax.set_ylim(-0.8, 0.8)
        ax.grid(True)
        ax.legend(fontsize=8)

        # Draw legend if at least one plotted element has label
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(fontsize=8, loc='upper left', framealpha=0.7)

figTitle = f"ALU_vs_phi_2D_{current_timestamp()}.png"
plt.savefig(figTitle, dpi=300)
plt.show()
print(f"ALU Binned Histogram saved as {figTitle}.")
