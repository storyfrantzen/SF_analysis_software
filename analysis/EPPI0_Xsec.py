# type: ignore
import argparse
import ROOT
import numpy as np
import math
from scipy.optimize import curve_fit

ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

# constants:
ELECTRON_CHARGE = 1.602e-19 # Coulombs
TARGET_LENGTH = 5.00 # cm
RHO_LH2 = 0.07 # g / cm3
N_A = 6.02214e23
ALPHA_EM = 1/137.035999084
PROTON_MASS = 0.9382720813  # proton mass [GeV]
PI = math.pi

# physics-motivated fit windows:
DEFAULT_FITS = {
    "m_gg":       {"fit_func": "gaus(0)+pol1(3)",  "range": (0.098, 0.17), "p0": None},      # pi0
    "m2_epX":     {"fit_func": "gaus(0)+pol1(3)",  "range": (-0.2, 0.2),   "p0": None},      # missing mass^2 epX
    "m2_epi0X":   {"fit_func": "gaus(0)+pol1(3)",  "range": (0.5, 1.5),    "p0": None},      # missing mass^2 epi0X
    "m_eggX":     {"fit_func": "gaus(0)+pol1(3)",  "range": (0.7, 1.4),    "p0": None},      # missing mass eggX
    "E_miss":     {"fit_func": "gaus(0)+pol1(3)",  "range": (-0.2, 0.4),   "p0": None},      # missing energy
    "pi0_thetaX": {"fit_func": "crystalball(0)",   "range": (0, 0.12),     "p0": None}       # pi0 thetaX
}

### ----------------- Helper Functions ----------------- ###
def coerce_scalar_column(arr):
    out = []
    for ev in arr:
        if ev is None:
            out.append(np.nan)
        elif isinstance(ev, (bytes, str, np.bytes_)):
            try: out.append(float(ev))
            except: out.append(np.nan)
        elif hasattr(ev, "__len__") and not isinstance(ev, (float,int)):
            try: out.append(float(ev[0]) if len(ev) > 0 else np.nan)
            except: out.append(np.nan)
        else:
            try: out.append(float(ev))
            except: out.append(np.nan)
    return np.array(out, dtype=float)

def get_adaptive_edges(tree, varname, n_bins, min_val, max_val, fine_bins=1000):
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

def get_bin_indices(values, edges):
    vals = np.asarray(values, dtype=float)
    indices = np.searchsorted(edges, vals, side='right') - 1
    valid = (indices >= 0) & (indices < len(edges)-1)
    indices = indices.astype(int)
    indices[~valid] = -1
    return indices

def fit_exclusive(vals, name=None, fit_func=None, fit_range=None, p0=None, nbins=None, fallback=False):
    """
    Fit a 1D distribution and return the mean and sigma.
    
    Parameters
    ----------
    vals : array-like
        Data points to fit.
    name : str
        Optional name to look up defaults from DEFAULT_FITS.
    fit_func : ROOT.TF1, optional
        Fit function (overrides DEFAULT_FITS if provided).
    fit_range : tuple, optional
        Fit range (overrides DEFAULT_FITS if provided).
    p0 : list/tuple, optional
        Initial parameter guesses.
    nbins : int, optional
        Number of bins for the histogram.
    fallback : bool
        If True, returns mean/std if fit fails.
    
    Returns
    -------
    mu, sigma : float
        Fitted mean and sigma.
    """
    vals = np.asarray(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    n = len(vals)
    if n == 0:
        return np.nan, np.nan

    if nbins is None:
        nbins = max(30, min(200, int(max(10, n // 2))))
    if n < 10:
        if fallback:
            mu, sigma = np.mean(vals), np.std(vals, ddof=1) if n > 1 else np.nan
            return mu, sigma
        return np.nan, np.nan

    # Histogram
    hist_name = f"h_fit_{np.random.randint(1_000_000_000)}"
    hist = ROOT.TH1D(hist_name, hist_name, nbins, np.min(vals), np.max(vals))
    for v in vals: hist.Fill(float(v))

    data_min, data_max = np.min(vals), np.max(vals)

    # Determine default fit settings from dict
    default_fit = DEFAULT_FITS.get(name, {})
    func_str = fit_func or default_fit.get("fit_func", "gaus(0)+pol1(3)")
    fit_range_default  = fit_range or default_fit.get("range", (data_min, data_max))
    p0_default = default_fit.get("p0", None)

    # Adjust fit range to actual data
    if fit_range is None:
        fit_min = max(fit_range_default[0], data_min)
        fit_max = min(fit_range_default[1], data_max)
        if fit_min >= fit_max:
            # fallback if no overlap
            if fallback:
                return np.mean(vals), np.std(vals, ddof=1)
            return np.nan, np.nan
        fit_range = (fit_min, fit_max)

    # Create ROOT TF1
    f = fit_func.Clone() if fit_func else ROOT.TF1(f"f_fit_{hist_name}", func_str, *fit_range)

    # Set initial parameters
    if p0 is None:
        if p0_default is not None:
            p0 = p0_default
        else:
            amp   = hist.GetMaximum()
            peak  = hist.GetBinCenter(hist.GetMaximumBin())
            width = hist.GetStdDev() if hist.GetStdDev() > 0 else 0.05
            if "crystalball" in func_str.lower():
                # amp, μ, σ, alpha (-1 * where the tail begins), n (tail power)
                p0 = (amp, peak, width, -1.5, 4)
            else:  # gaussian+pol1 default
                p0 = (amp, peak, width, 0, 0)

    for i, val in enumerate(p0):
        f.SetParameter(i, val)

    # Perform fit quietly
    hist.Fit(f, "RQ0")

    # Extract mu and sigma based on function type
    try:
        title = f.GetTitle().lower()
        if "crystalball" in title:
            mu, sigma = f.GetParameter(1), abs(f.GetParameter(2))
        else:  # gaussian-style
            mu, sigma = f.GetParameter(1), abs(f.GetParameter(2))
    except Exception:
        if fallback:
            mu, sigma = np.mean(vals), np.std(vals, ddof=1)
        else:
            mu, sigma = np.nan, np.nan

    hist.Delete()
    return mu, sigma

def load_gemc(gemc_file, cols):
    """Load GEMC ROOT file into numpy dict."""
    f_gemc = ROOT.TFile(gemc_file, "READ")
    tree = f_gemc.Get("Events")
    df_gemc = ROOT.RDataFrame(tree)

    df_numpy = df_gemc.AsNumpy(columns=cols)
    for k in cols:
        df_numpy[k] = coerce_scalar_column(df_numpy[k])
    return df_numpy

def compute_acceptance(df_numpy_gemc, Q2_edges, Xb_edges, t_edges, phi_edges, shape):
    """Return acceptance[Q2, Xb, t, phi] from GEMC numpy arrays."""
    gen_counts = np.zeros(shape, dtype=float)
    rec_counts = np.zeros(shape, dtype=float)

    iQ2_gen  = get_bin_indices(df_numpy_gemc["gen_dis.Q2"], Q2_edges)
    iXb_gen  = get_bin_indices(df_numpy_gemc["gen_dis.Xb"], Xb_edges)
    it_gen   = get_bin_indices(df_numpy_gemc["gen_eppi0.t"],  t_edges)
    iphi_gen = get_bin_indices((df_numpy_gemc["gen_eppi0.trentoPhi"]+2*np.pi)%(2*np.pi)*180/np.pi, phi_edges)

    iQ2_rec  = get_bin_indices(df_numpy_gemc["dis.Q2"], Q2_edges)
    iXb_rec  = get_bin_indices(df_numpy_gemc["dis.Xb"], Xb_edges)
    it_rec   = get_bin_indices(df_numpy_gemc["eppi0.t"],  t_edges)
    iphi_rec = get_bin_indices((df_numpy_gemc["eppi0.trentoPhi"]+2*np.pi)%(2*np.pi)*180/np.pi, phi_edges)

    n_events = len(df_numpy_gemc["gen_dis.Q2"])
    for i in range(n_events):
        # generated
        if iQ2_gen[i] >= 0 and iXb_gen[i] >= 0 and it_gen[i] >= 0 and iphi_gen[i] >= 0:
            gen_counts[iQ2_gen[i], iXb_gen[i], it_gen[i], iphi_gen[i]] += 1
        # reconstructed
        if iQ2_rec[i] >= 0 and iXb_rec[i] >= 0 and it_rec[i] >= 0 and iphi_rec[i] >= 0:
            rec_counts[iQ2_rec[i], iXb_rec[i], it_rec[i], iphi_rec[i]] += 1

    acceptance = np.divide(rec_counts, gen_counts, out=np.zeros_like(rec_counts), where=gen_counts>0)
    return acceptance

def save_acceptance_maps(acceptance, Q2_edges, Xb_edges, t_edges, phi_edges, filename="acceptance_maps.root"):
    """Write acceptance histograms to ROOT file for inspection."""
    out_file = ROOT.TFile(filename,"RECREATE")
    for q2 in range(acceptance.shape[0]):
        for xb in range(acceptance.shape[1]):
            for t in range(acceptance.shape[2]):
                h_acc = ROOT.TH1D(
                    f"h_acc_q{q2}_xb{xb}_t{t}",
                    f"Acceptance Q2[{Q2_edges[q2]:.2f}-{Q2_edges[q2+1]:.2f}], "
                    f"Xb[{Xb_edges[xb]:.2f}-{Xb_edges[xb+1]:.2f}], "
                    f"-t[{t_edges[t]:.2f}-{t_edges[t+1]:.2f}]",
                    len(phi_edges)-1, phi_edges
                )
                for phi in range(len(phi_edges)-1):
                    h_acc.SetBinContent(phi+1, acceptance[q2, xb, t, phi])
                h_acc.Write()
    out_file.Close()

def epsilon(Q2, y, E):
    numerator = 1 - y - Q2/(4*E**2)
    denom = 1 - y + y**2/2 + Q2/(4*E**2)
    if denom <= 0:
        return 0.0
    return numerator / denom

def get_y(Q2, xb, E):
    return Q2 / (2.0 * PROTON_MASS * xb * E)

def gamma_flux(Q2, xb, E):
    if Q2 <= 0 or xb <= 0:
        return 0.0

    # compute epsilon
    y = get_y(xb, Q2, E)
    eps = epsilon(y, Q2, E)
    # protect against eps >= 1 (would blow up 1/(1-eps))
    if eps >= 0.999999:
        return 0.0
    
    flux = ALPHA_EM / (8.0 * PI) * Q2 / (PROTON_MASS * E)**2 * (1 - xb) / xb / xb / xb * 1 / (1 - eps)

    # safety: require flux finite and non-negative
    if not np.isfinite(flux) or flux <= 0:
        return 0.0
    return flux

def luminosity(q_beam, l_target=TARGET_LENGTH, rho_target=RHO_LH2):
    # returns luminosity in inverse femtobarns, provided input beam charge in C
    return N_A * l_target * rho_target * q_beam / ELECTRON_CHARGE * 1e-39  # fb^-1

### ----------------- ARGS ----------------- ###
parser = argparse.ArgumentParser(description="Compute XSection from input ROOT file.")
parser.add_argument("input_file", type=str)
parser.add_argument("E", type=float, help="Beam energy")
parser.add_argument("--gemc", type=str, help="Optional GEMC file for acceptance corrections")
parser.add_argument("-a","--adaptive", action="store_true")
parser.add_argument("-c6", "--clas6", action="store_true")
parser.add_argument("-c12", "--clas12", action="store_true")
parser.add_argument("-m", "--manual", action="store_true")
parser.add_argument("--global_cuts", action="store_true", help="Derive exclusivity cuts globally instead of bin-by-bin")
parser.add_argument("-v", "--verbose", action="store_true", help="Verbose, diagnostics for all kinematic bins.")
args = parser.parse_args()

### --------------- Read eppi0REC ROOT file --------------- ###
f = ROOT.TFile(args.input_file, "READ")
summary_tree = f.Get("Summary")
summary_tree.GetEntry(0)
events_tree = f.Get("Events")
df = ROOT.RDataFrame(events_tree)
cols = ["dis.Q2","dis.Xb","eppi0.t","eppi0.trentoPhi","p.det","eppi0.pi0_thetaX","eppi0.m2_epX","eppi0.m_eggX","eppi0.E_miss","eppi0.m_gg"]
df_numpy = df.AsNumpy(columns=cols)
for k in cols: df_numpy[k] = coerce_scalar_column(df_numpy[k])
ex_vars = ["m2_epX", "m_eggX", "E_miss", "pi0_thetaX", "m_gg"]
exclusive_vars = ["m2_epX", "m_eggX", "E_miss", "pi0_thetaX"] # not including m_gg

### --------------- Initialize binning scheme --------------- ###
if args.adaptive:
    Q2_edges = get_adaptive_edges(events_tree, "dis.Q2", 6, 1.0, 6.5)
    Xb_edges = get_adaptive_edges(events_tree, "dis.Xb", 6, 0.1, 0.7)
    t_edges  = get_adaptive_edges(events_tree, "eppi0.t", 6, 0.0, 2.0)
    print(f"Adaptive binning yielded {len(Q2_edges)} Q2 bins, {len(Xb_edges)} Xb bins, {len(t_edges)} t bins.")
    phi_edges = np.linspace(0, 360, 21)
elif args.clas6:
    Q2_edges = np.array([1.0,1.5,2,2.5,3,3.5,4,4.6])
    Xb_edges = np.array([0.1,0.15,0.2,0.25,0.3,0.38,0.48,0.58])
    t_edges  = np.array([0.09,0.15,0.2,0.3,0.4,0.6,1.0,1.5,2.0])
    phi_edges = np.linspace(0, 360, 21)
elif args.clas12:
    Q2_edges = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4, 4.6, 5.5, 7.0, 10.5]) # guessed final three edges from figures
    Xb_edges = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.38, 0.48, 0.58, 0.9])
    t_edges  = np.array([0.09, 0.15, 0.2, 0.3, 0.4, 0.6, 1.0, 1.5, 2.0])
    phi_edges = np.linspace(0, 360, 11)
elif args.manual:
    Q2_edges = np.array([1.0, 1.5, 2, 2.5, 3, 3.5, 4, 6.5])
    #Q2_edges = np.array([1.0, 7.5])
    Xb_edges = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
    t_edges  = np.array([0.0, 2.0])
    phi_edges = np.linspace(0, 360, 11)
else:
    raise RuntimeError("Must specify binning scheme, e.g., --adaptive or --clas6 or --manual")

iQ2   = get_bin_indices(df_numpy["dis.Q2"], Q2_edges)
iXb   = get_bin_indices(df_numpy["dis.Xb"], Xb_edges)
it    = get_bin_indices(df_numpy["eppi0.t"], t_edges)
iphi  = get_bin_indices((df_numpy["eppi0.trentoPhi"] + 2*np.pi) % (2*np.pi) * 180.0/np.pi, phi_edges)
pIdx  = np.where(np.isclose(df_numpy["p.det"], 1.0), 0, 1).astype(int)

in_range = (iQ2>=0) & (iXb>=0) & (it>=0) & (iphi>=0)
print(f"Total events in kinematic range: {np.sum(in_range)}")

bins_map = (len(Q2_edges)-1, len(Xb_edges)-1, len(t_edges)-1, len(phi_edges)-1, 2)

### ------ Sideband subtraction parameters:
n_sigma_signal = 3
n_sigma_sb_min = 3
n_sigma_sb_max = 5

### -------- Global ex. cuts (per proton topology) --------- ###
global_cuts = {}
if args.global_cuts:
    print(">>> Deriving global exclusivity cuts (for separate proton topologies)...")
    for p in range(2):  # pIdx = 0 (FD), 1 (CD)
        surv_mask = (in_range) & (pIdx == p)
        for var in exclusive_vars:
            vals = df_numpy[f"eppi0.{var}"][surv_mask]
            if len(vals) == 0:
                continue
            mu, sigma = fit_exclusive(vals, name=var)
            global_cuts[(var, p)] = (mu, sigma)
            if args.verbose:
                print(f"Global {var}, pIdx={p}: mu={mu:.5f}, sigma={sigma:.5f}")

             # apply cut immediately so next var sees reduced distribution
            surv_mask[surv_mask] &= np.abs(vals - mu) <= n_sigma_signal*sigma
            if args.verbose:
                print(f"Var={var}, pIdx={p}: mu={mu:.5f}, sigma={sigma:.5f}, survived {np.sum(surv_mask)}")

### ------ create 4D arrays + boolean mask for final survivors:
side_sub_yield4D = np.zeros(bins_map[:4], dtype=float) 
errs4D = np.zeros(bins_map[:4], dtype=float)
xsec4D = np.zeros(bins_map[:4], dtype=float) 
surv_mask_all = np.zeros(len(df_numpy["dis.Q2"]), dtype=bool)

BEAM_Q = summary_tree.TotalCharge * 1e-9 # Coulombs
LUM_INT = luminosity(BEAM_Q) # fb^-1
BR = 0.988

### ----------------- Acceptance from GEMC: ----------------- ###
acceptance = None
if args.gemc:
    cols_gemc = ["gen_dis.Q2","gen_dis.Xb","gen_eppi0.t","gen_eppi0.trentoPhi", "dis.Q2","dis.Xb","eppi0.t","eppi0.trentoPhi"]
    df_numpy_gemc = load_gemc(args.gemc, cols_gemc)
    acceptance = compute_acceptance(df_numpy_gemc, Q2_edges, Xb_edges, t_edges, phi_edges, xsec4D.shape)
    save_acceptance_maps(acceptance, Q2_edges, Xb_edges, t_edges, phi_edges)

### ----------- Create tree to track raw counts ----------- ###
out_file = ROOT.TFile("survival_summary.root", "RECREATE")
tree_raw = ROOT.TTree("raw_counts", "Raw event counts per bin")
tree_raw.SetAutoFlush(10000)  # flush every 10k entries

q2b   =  np.zeros(1, dtype=int)
xbb   = np.zeros(1, dtype=int)
tb    = np.zeros(1, dtype=int)
phib  = np.zeros(1, dtype=int)
pIdxb = np.zeros(1, dtype=int)
count = np.zeros(1, dtype=int)

tree_raw.Branch("q2_bin", q2b, "q2_bin/I")
tree_raw.Branch("xb_bin", xbb, "xb_bin/I")
tree_raw.Branch("t_bin", tb, "t_bin/I")
tree_raw.Branch("phi_bin", phib, "phi_bin/I")
tree_raw.Branch("count", count, "count/I")

### ---------------------------------- XSECTION Workflow ---------------------------------- ###
for q2 in range(bins_map[0]):
    dq2 = Q2_edges[q2+1] - Q2_edges[q2]
    for xb in range(bins_map[1]):
        dxb = Xb_edges[xb+1] - Xb_edges[xb]
        for t in range(bins_map[2]):
            dt = t_edges[t+1] - t_edges[t]
            for phi in range(bins_map[3]):
                dphi = phi_edges[phi+1] - phi_edges[phi]
        
                mask_bin = (in_range) & (iQ2==q2) & (iXb==xb) & (it==t) & (iphi==phi)
                num_events = np.sum(mask_bin)
                if num_events > 0:
                    q2b[0], xbb[0], tb[0], phib[0], count[0] = q2, xb, t, phi, num_events
                    tree_raw.Fill()

                bin_volume = dq2 * dxb * dt * dphi
                bin_variance = 0
                bin_yield = 0

                # loop over 2 possible proton topologies:
                for p in range(2):
                    mask_bin_p = mask_bin & (pIdx==p)
                    num_events = np.sum(mask_bin_p)
                    
                    if args.verbose and num_events > 0:
                        print(f"Bin [Q2: {q2}, Xb: {xb}, t: {t}, phi: {phi}, pIdx: {p}] += {num_events} events.")
                    if not np.any(mask_bin_p):
                        continue
                    
                    surv_mask = mask_bin_p.copy()
                    for var in exclusive_vars:
                        if args.verbose and np.sum(surv_mask) > 0:
                            print(f"Initiating {var} cuts with {np.sum(surv_mask)} events.")
                        if args.global_cuts:
                            if (var, p) not in global_cuts:
                                continue  # no global fit possible for this topology
                            mu, sigma = global_cuts[(var, p)]
                        else:
                            vals = df_numpy[f"eppi0.{var}"][surv_mask]
                            mu, sigma = fit_exclusive(vals, name=var, nbins=200)

                        vals = df_numpy[f"eppi0.{var}"][surv_mask]
                        surv_mask[surv_mask] &= np.abs(vals - mu) <= n_sigma_signal*sigma

                        if args.verbose:
                            print(f"{var} cut applied with mu={mu:.5f}, sigma={sigma:.5f}, "
                                  f"survived {np.sum(surv_mask)} events")
                        
                    
                    ### ------------ Sideband subtraction (m_gg) ------------ ###
                    mgg_vals = df_numpy["eppi0.m_gg"][surv_mask]
                    mu_mgg, sigma_mgg = fit_exclusive(mgg_vals, name="m_gg", nbins=200)
                    sig_mask = np.abs(mgg_vals - mu_mgg) < n_sigma_signal*sigma_mgg
                    sb_mask  = ((mu_mgg - n_sigma_sb_max*sigma_mgg <= mgg_vals) & 
                                (mgg_vals < mu_mgg - n_sigma_sb_min*sigma_mgg)) | \
                               ((mu_mgg + n_sigma_sb_min*sigma_mgg < mgg_vals) & 
                                (mgg_vals <= mu_mgg + n_sigma_sb_max*sigma_mgg))

                    n_sig = np.sum(sig_mask)
                    n_sb  = np.sum(sb_mask)
                    w_sig = 2 * n_sigma_signal * sigma_mgg
                    w_sb  = 2 * (n_sigma_sb_max - n_sigma_sb_min) * sigma_mgg
                    alpha = w_sig / w_sb
                    if w_sb > 0:
                        n_bkg = n_sb * alpha
                        if args.verbose and n_sig > 0:
                            print(f"Bin [Q2: {q2}, Xb: {xb}, t: {t}, phi: {phi}] contains {n_sig} signal events and {n_bkg} background events.")
                        bin_variance += n_sig + (alpha**2) * n_sb
                        bin_yield += max(n_sig - n_bkg, 0.0)

                    # Mark surviving events
                    surv_mask_all[surv_mask] = True
                
                side_sub_yield4D[q2, xb, t, phi] = bin_yield
                errs4D[q2, xb, t, phi] = np.sqrt(bin_variance) / bin_volume / LUM_INT / BR
                xsec4D[q2, xb, t, phi] = bin_yield / bin_volume / LUM_INT / BR

if acceptance is not None:
    xsec4D = np.divide(xsec4D, acceptance, out=np.zeros_like(xsec4D), where=acceptance>0)
    errs4D = np.divide(errs4D, acceptance, out=np.zeros_like(errs4D), where=acceptance>0)

surviving_indices = np.flatnonzero(surv_mask_all)
print(f"Total surviving events: {len(surviving_indices)}")
print(f"Total yield after 3σ cuts + sideband subtraction: {side_sub_yield4D.sum():.3f}")

### ----------------- Fill survival ROOT histograms ----------------- ###
hist_dict = {}
graph_list = []
fit_list = []

nbins_default = 100

for var in ex_vars:
    vals = df_numpy[f"eppi0.{var}"][surviving_indices]
    if len(vals) == 0: 
        continue
    h = ROOT.TH1D(f"h_{var}_all", f"h_{var}_all", nbins_default, np.min(vals), np.max(vals))
    for v in vals: h.Fill(float(v))
    h.Write()
    hist_dict[var] = h

h_Q2_Xb = ROOT.TH2D("h_Q2_Xb", "Q^{2} vs x_{B} Coverage; x_{B}; Q^{2} [GeV^{2}]; Events",
                    len(Xb_edges)-1, Xb_edges, len(Q2_edges)-1, Q2_edges)
for q2, xb in zip(df_numpy["dis.Q2"][surviving_indices], df_numpy["dis.Xb"][surviving_indices]): h_Q2_Xb.Fill(xb, q2)
h_Q2_Xb.Write()

tree_raw.AutoSave("FlushBaskets")
out_file.Close()
print("Saved exclusive histograms to exclusive_summary.root")

### ----------------- Project 4D histogram to phi ----------------- ###
out_file = ROOT.TFile("phi_xsec.root", "RECREATE")
Ebeam = args.E

phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
phi_widths = 0.5 * (phi_edges[1:] - phi_edges[:-1])

for q2 in range(bins_map[0]):
    for xb in range(bins_map[1]):
        for t in range(bins_map[2]):
            phi_yields = xsec4D[q2, xb, t, :]
            phi_errors = errs4D[q2, xb, t, :]

            if np.sum(phi_yields) == 0:
                continue

            # surviving event mask for this 4D bin
            mask_bin = (iQ2 == q2) & (iXb == xb) & (it == t) & surv_mask_all
            if not np.any(mask_bin):
                continue

            # compute bin centers for Q2 and Xb (arithmetic mean of surviving events)
            Q2_center = np.mean(df_numpy["dis.Q2"][mask_bin])
            Xb_center = np.mean(df_numpy["dis.Xb"][mask_bin])

            # compute virtual photon flux
            gamma = gamma_flux(Q2_center, Xb_center, Ebeam)
            if gamma <= 0:
                continue

            # divide by flux to get reduced cross section
            xsec_red = phi_yields / gamma
            xsec_red_err = phi_errors / gamma

            n_points = len(phi_centers)
            gr_name = f"gr_phi_q{q2}_xb{xb}_t{t}"
            gr_title = (
                f"Q2={Q2_edges[q2]:.2f}-{Q2_edges[q2+1]:.2f}, "
                f"Xb={Xb_edges[xb]:.2f}-{Xb_edges[xb+1]:.2f}, "
                f"-t={t_edges[t]:.2f}-{t_edges[t+1]:.2f}; "
                f"#phi [deg]; #pi^0 Reduced Cross Section [fb/GeV^2]"
            )

            gr = ROOT.TGraphErrors(
                n_points,
                phi_centers.astype(np.float64),
                xsec_red.astype(np.float64),
                phi_widths.astype(np.float64),
                xsec_red_err.astype(np.float64)
            )

            # fit: A + B cos(phi) + C cos(2 phi)
            fit_func = ROOT.TF1("fit_phi", "[0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad())", 0, 360)
            fit_func.SetParameters(np.mean(xsec_red), 0.1*np.mean(xsec_red), 0.1*np.mean(xsec_red))
            gr.Fit(fit_func, "QR")

            gr.SetName(gr_name)
            gr.SetTitle(gr_title)
            gr.Draw("AP")
            fit_func.Draw("same")
            gr.GetXaxis().SetRangeUser(0, 360)
            gr.Write()

            graph_list.append(gr)      # keep reference
            fit_list.append(fit_func)  # keep reference

out_file.Close()
print("Saved phi-projected reduced cross sections to phi_xsec.root")