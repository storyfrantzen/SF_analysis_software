# type: ignore
import argparse
import ROOT
import numpy as np

ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

### ----------------- Helper Functions ----------------- ###
def get_adaptive_edges(tree, varname, n_bins, min_val, max_val, fine_bins=500):
    hist_name = f"h_{varname}"
    hist = ROOT.TH1D(hist_name, hist_name, fine_bins, min_val, max_val)
    tree.Draw(f"{varname} >> {hist_name}", "", "goff")
    total_entries = hist.GetEntries()
    if total_entries == 0:
        hist.Delete()
        return np.linspace(min_val, max_val, n_bins+1)
    step = total_entries / n_bins
    edges = [hist.GetXaxis().GetXmin()]
    cumulative, target = 0, step
    for i in range(1, hist.GetNbinsX() + 1):
        cumulative += hist.GetBinContent(i)
        if cumulative >= target:
            edge = hist.GetBinLowEdge(i + 1)
            if edge > edges[-1]:
                edges.append(edge)
                target += step
    if edges[-1] < hist.GetXaxis().GetXmax():
        edges.append(hist.GetXaxis().GetXmax())
    hist.Delete()
    return np.array(edges)

def fit_exclusive(vals, fit_func=None, fit_range=None, p0=None, nbins=None):
    vals = np.asarray(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    n = len(vals)
    if n == 0:
        return np.nan, np.nan
    if nbins is None:
        nbins = max(30, min(200, int(max(10, n // 2))))
    if n < 10:
        mu, sigma = np.mean(vals), np.std(vals, ddof=1) if n > 1 else np.nan
        return mu, sigma
    hist_name = f"h_fit_{np.random.randint(1_000_000_000)}"
    hist = ROOT.TH1D(hist_name, hist_name, nbins, np.min(vals), np.max(vals))
    for v in vals: hist.Fill(float(v))
    f = fit_func.Clone() if fit_func else ROOT.TF1("f_fit_temp","gaus(0)+pol1(3)", np.min(vals), np.max(vals))
    if fit_range: f.SetRange(*fit_range)
    if p0:
        for i, val in enumerate(p0):
            f.SetParameter(i, val)
    hist.Fit(f, "RQ0")
    try:
        mu, sigma = f.GetParameter(1), abs(f.GetParameter(2))
    except Exception:
        mu, sigma = np.mean(vals), np.std(vals, ddof=1)
    hist.Delete()
    return mu, sigma

def coerce_scalar_column(arr):
    """
    Convert AsNumpy column into a 1D float array of length nEvents.
    If arr elements are arrays/lists, take first element; if empty -> nan.
    """
    out = []
    for ev in arr:
        # handle ROOT vector types (they come through as numpy.ndarray or list-like)
        if ev is None:
            out.append(np.nan)
        else:
            # bytes/str special-case: try to cast
            if isinstance(ev, (bytes, str, np.bytes_)):
                try:
                    out.append(float(ev))
                except Exception:
                    out.append(np.nan)
            # ndarray or list-like
            elif hasattr(ev, "__len__") and not isinstance(ev, (float, int)):
                try:
                    out.append(float(ev[0]) if len(ev) > 0 else np.nan)
                except Exception:
                    out.append(np.nan)
            else:
                try:
                    out.append(float(ev))
                except Exception:
                    out.append(np.nan)
    return np.array(out, dtype=float)

def get_bin_indices(values, edges):
    vals = np.asarray(values, dtype=float)
    indices = np.searchsorted(edges, vals, side='right') - 1
    valid = (indices >= 0) & (indices < len(edges)-1)
    indices = indices.astype(int)
    indices[~valid] = -1
    return indices

### ----------------- Main Workflow ----------------- ###
parser = argparse.ArgumentParser(description="Compute XSection from input ROOT file.")
parser.add_argument("input_file", type=str)
parser.add_argument("--adaptive", action="store_true")
parser.add_argument("--clas6", action="store_true")
parser.add_argument("--plot", action="store_true", help="Plot histogram slices")
args = parser.parse_args()

# Load data
f = ROOT.TFile(args.input_file, "READ")
tree = f.Get("Events")
df = ROOT.RDataFrame(tree)

cols_needed = ["dis.Q2","dis.Xb","eppi0.t","eppi0.trentoPhi","p.det",
               "eppi0.pi0_thetaX","eppi0.m2_epX","eppi0.m_eggX","eppi0.E_miss","eppi0.m_gg"]

df_numpy = df.AsNumpy(columns=cols_needed)

# Coerce to scalars safely (prevents casting issues if input tree not setup correctly)
for k in cols_needed:
    df_numpy[k] = coerce_scalar_column(df_numpy[k])

# Bin edges
if args.adaptive:
    Q2_edges = get_adaptive_edges(tree, "dis.Q2", 7, 1.0, 6.5)
    Xb_edges = get_adaptive_edges(tree, "dis.Xb", 7, 0.1, 0.7)
    t_edges  = get_adaptive_edges(tree, "eppi0.t", 8, 0.0, 4.0)
elif args.clas6:
    Q2_edges = np.array([1.0,1.5,2,2.5,3,3.5,4,4.6])
    Xb_edges = np.array([0.1,0.15,0.2,0.25,0.3,0.38,0.48,0.58])
    t_edges  = np.array([0.09,0.15,0.2,0.3,0.4,0.6,1.0,1.5,2.0])
else:
    raise RuntimeError("Must specify --adaptive or --clas6")

phi_edges = np.linspace(0, 360, 21)

# Compute 4D bin indices
iQ2   = get_bin_indices(df_numpy["dis.Q2"], Q2_edges)
iXb   = get_bin_indices(df_numpy["dis.Xb"], Xb_edges)
it    = get_bin_indices(df_numpy["eppi0.t"], t_edges)
trento = np.nan_to_num(df_numpy["eppi0.trentoPhi"], nan=0.0)
iphi  = get_bin_indices((trento + 2*np.pi) % (2*np.pi) * 180.0/np.pi, phi_edges)
pdet = df_numpy["p.det"]
pIdx  = np.where(np.isclose(pdet, 1.0), 0, 1).astype(int)

# --- valid mask ---
in_range = (iQ2>=0) & (iXb>=0) & (it>=0) & (iphi>=0)
n_events = np.sum(in_range)
print(f"Total events in kinematic range: {n_events}")

bins_map = (len(Q2_edges)-1, len(Xb_edges)-1, len(t_edges)-1, len(phi_edges)-1, 2)
ex_vars = ["pi0_thetaX","m2_epX","m_eggX","E_miss","m_gg"]

# Initialize mu/sigma with NaN so empty bins are explicit
mu_dict = {v: np.full(bins_map, np.nan, dtype=float) for v in ex_vars}
sigma_dict = {v: np.full(bins_map, np.nan, dtype=float) for v in ex_vars}

# Fit mu/sigma for each bin (only using events that fall in that bin)
for q in range(bins_map[0]):
    for x in range(bins_map[1]):
        for t in range(bins_map[2]):
            for phi in range(bins_map[3]):
                for p in range(2):
                    mask = (iQ2==q) & (iXb==x) & (it==t) & (iphi==phi) & (pIdx==p) & in_range
                    if not np.any(mask): 
                        continue
                    for var in ex_vars:
                        vals = df_numpy[f"eppi0.{var}"][mask]
                        mu, sigma = fit_exclusive(vals)
                        mu_dict[var][q,x,t,phi,p] = mu
                        sigma_dict[var][q,x,t,phi,p] = sigma

# ----------------- Background subtraction and final histogram (with 3σ cuts) -----------------
n_sigma_signal = 3
n_sigma_sb_min = 3
n_sigma_sb_max = 6

hist4d_bkgsub = np.zeros(bins_map[:4], dtype=float)  # background-subtracted counts BEFORE 3σ on other vars
hist4d_final  = np.zeros(bins_map[:4], dtype=float)  # background-subtracted counts AFTER 3σ

# --- Diagnostic mode for per-bin statistics ---
diagnostic = False  # set False to skip

# First: compute hist4d_before (sideband subtraction using all events in bin)
for q2 in range(bins_map[0]):
    for xb in range(bins_map[1]):
        for t in range(bins_map[2]):
            for phi in range(bins_map[3]):
                for p in range(2):
                    mask = (iQ2==q2) & (iXb==xb) & (it==t) & (iphi==phi) & (pIdx==p) & in_range
                    n_events = np.sum(mask)
                    if n_events == 0:
                        if diagnostic:
                            print(f"Empty bin: Q2={q2}, Xb={xb}, t={t}, phi={phi}, p={p}")
                        continue

                    mgg = df_numpy["eppi0.m_gg"][mask]
                    mu = mu_dict["m_gg"][q2,xb,t,phi,p]
                    sigma = sigma_dict["m_gg"][q2,xb,t,phi,p]

                    if not np.isfinite(mu) or not np.isfinite(sigma) or sigma <= 0:
                        continue

                    sig_mask = np.abs(mgg - mu) < n_sigma_signal*sigma
                    sb_mask = ((mu - n_sigma_sb_max*sigma <= mgg) & (mgg < mu - n_sigma_sb_min*sigma)) | \
                              ((mu + n_sigma_sb_min*sigma < mgg) & (mgg <= mu + n_sigma_sb_max*sigma))
                    n_sig = np.sum(sig_mask)
                    n_sb  = np.sum(sb_mask)
                    w_sig = 2*n_sigma_signal*sigma
                    w_sb  = 2*(n_sigma_sb_max - n_sigma_sb_min)*sigma

                    if w_sb <= 0:
                        continue
                    n_bkg = n_sb * (w_sig / w_sb)
                    hist4d_bkgsub[q2,xb,t,phi] += max(n_sig - n_bkg, 0.0)

                    if diagnostic:
                        print(f"Bin Q2={q2}, Xb={xb}, t={t}, phi={phi}, p={p}, n_events={n_events}, n_sig={n_sig}, n_bkg={n_bkg}, contrib={max(n_sig-n_bkg,0)}")

print(f"Sum of hist4d_before (no 3σ on other vars): {hist4d_bkgsub.sum():.3f}")

# Build final histogram by applying the other 3σ cuts *per event* inside each bin,
# then doing sideband subtraction only on the surviving events in that bin.
for q in range(bins_map[0]):
    for x in range(bins_map[1]):
        for t in range(bins_map[2]):
            for phi in range(bins_map[3]):
                for p in range(2):
                    mask_bin = (iQ2==q) & (iXb==x) & (it==t) & (iphi==phi) & (pIdx==p) & in_range
                    if not np.any(mask_bin):
                        continue
                    # collect mu & sigma for the four other vars
                    mus = {}
                    sigmas = {}
                    skip_bin = False
                    for var in ["pi0_thetaX","m2_epX","m_eggX","E_miss"]:
                        mu = mu_dict[var][q,x,t,phi,p]
                        sigma = sigma_dict[var][q,x,t,phi,p]
                        if not np.isfinite(mu) or not np.isfinite(sigma) or sigma <= 0:
                            skip_bin = True
                            break
                        mus[var] = mu
                        sigmas[var] = sigma
                    if skip_bin:
                        continue
                    events_idx = np.where(mask_bin)[0]
                    # apply 3σ cuts for all listed vars to keep survivors
                    keep_mask = np.ones(len(events_idx), dtype=bool)
                    for var in ["pi0_thetaX","m2_epX","m_eggX","E_miss"]:
                        arr = df_numpy[f"eppi0.{var}"][events_idx]
                        keep_mask &= (np.abs(arr - mus[var]) <= 3.0*sigmas[var])
                    if not np.any(keep_mask):
                        continue
                    # now perform sideband subtraction on surviving events only
                    mgg_surv = df_numpy["eppi0.m_gg"][events_idx][keep_mask]
                    mu_mgg = mu_dict["m_gg"][q,x,t,phi,p]
                    sigma_mgg = sigma_dict["m_gg"][q,x,t,phi,p]
                    if not np.isfinite(mu_mgg) or not np.isfinite(sigma_mgg) or sigma_mgg <= 0:
                        continue
                    sig_mask = np.abs(mgg_surv - mu_mgg) < n_sigma_signal*sigma_mgg
                    sb_mask = ((mu_mgg - n_sigma_sb_max*sigma_mgg <= mgg_surv) & (mgg_surv < mu_mgg - n_sigma_sb_min*sigma_mgg)) | \
                              ((mu_mgg + n_sigma_sb_min*sigma_mgg < mgg_surv) & (mgg_surv <= mu_mgg + n_sigma_sb_max*sigma_mgg))
                    n_sig = np.sum(sig_mask)
                    n_sb  = np.sum(sb_mask)
                    w_sig = 2*n_sigma_signal*sigma_mgg
                    w_sb  = 2*(n_sigma_sb_max - n_sigma_sb_min)*sigma_mgg
                    if w_sb <= 0:
                        continue
                    n_bkg = n_sb * (w_sig / w_sb)
                    hist4d_final[q,x,t,phi] += max(n_sig - n_bkg, 0.0)

print(f"Sum after applying 3σ cuts (final): {hist4d_final.sum():.3f}")

# --- Visualization (optional) ---
if args.plot:
    import matplotlib.pyplot as plt
    # Example: project 4D histogram into Q2-Xb plane (summing over t & phi)
    proj_Q2Xb = hist4d_final.sum(axis=(2,3))
    plt.figure(figsize=(8,6))
    plt.imshow(proj_Q2Xb, origin='lower', aspect='auto',
               extent=[Xb_edges[0], Xb_edges[-1], Q2_edges[0], Q2_edges[-1]],
               cmap='viridis')
    plt.xlabel("Xb")
    plt.ylabel("Q2 [GeV^2]")
    plt.colorbar(label="Counts")
    plt.title("Projected 4D histogram (Q2-Xb)")
    plt.savefig("Xsec.png")
    plt.show()