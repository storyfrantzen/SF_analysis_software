# type: ignore
import sys
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
def is_valid_rec(df_numpy):
    mask = np.ones(len(df_numpy["dis.Q2"]), dtype=bool)
    for k, arr in df_numpy.items():
        if k.startswith("eppi0.") or k.startswith("dis."):
            if arr.dtype.kind in ("i", "u"):  # integers
                mask &= (arr != -999)
            else:  # floats
                mask &= ~np.isnan(arr)
    return mask

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
        # print("Fewer than 10 events, can't fit here...")
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

def compute_acceptance(df_numpy_gemc, Q2_edges, Xb_edges, t_edges, phi_edges, shape, valid_mask=None):
    """Return acceptance[Q2, Xb, t, phi] from GEMC numpy arrays."""
    gen_counts = np.zeros(shape, dtype=float)
    rec_counts = np.zeros(shape, dtype=float)

    # indices for generated
    iQ2_gen  = get_bin_indices(df_numpy_gemc["gen_dis.Q2"], Q2_edges)
    iXb_gen  = get_bin_indices(df_numpy_gemc["gen_dis.Xb"], Xb_edges)
    it_gen   = get_bin_indices(df_numpy_gemc["gen_eppi0.t"],  t_edges)
    iphi_gen = get_bin_indices((df_numpy_gemc["gen_eppi0.trentoPhi"]+2*np.pi)%(2*np.pi)*180/np.pi, phi_edges)

    # indices for reconstructed
    iQ2_rec  = get_bin_indices(df_numpy_gemc["dis.Q2"], Q2_edges)
    iXb_rec  = get_bin_indices(df_numpy_gemc["dis.Xb"], Xb_edges)
    it_rec   = get_bin_indices(df_numpy_gemc["eppi0.t"],  t_edges)
    iphi_rec = get_bin_indices((df_numpy_gemc["eppi0.trentoPhi"]+2*np.pi)%(2*np.pi)*180/np.pi, phi_edges)

    n_events = len(df_numpy_gemc["gen_dis.Q2"])
    for i in range(n_events):
        # --- generated always counts ---
        if iQ2_gen[i] >= 0 and iXb_gen[i] >= 0 and it_gen[i] >= 0 and iphi_gen[i] >= 0:
            gen_counts[iQ2_gen[i], iXb_gen[i], it_gen[i], iphi_gen[i]] += 1

        # --- reconstructed only if valid ---
        if valid_mask is None or valid_mask[i]:
            if iQ2_rec[i] >= 0 and iXb_rec[i] >= 0 and it_rec[i] >= 0 and iphi_rec[i] >= 0:
                rec_counts[iQ2_rec[i], iXb_rec[i], it_rec[i], iphi_rec[i]] += 1

    acceptance = np.divide(rec_counts, gen_counts, out=np.zeros_like(rec_counts), where=gen_counts>0)
    return acceptance, gen_counts, rec_counts

def save_acceptance_maps(acceptance, gen_counts, rec_counts,
                         Q2_edges, Xb_edges, t_edges, phi_edges,
                         filename="acceptance_maps.root"):
    """Write acceptance histograms and gen/rec overlays to ROOT file with proper canvases."""
    out_file = ROOT.TFile(filename, "RECREATE")

    # Global histogram of acceptance correction factors
    h_global = ROOT.TH1D("h_acceptance_all",
                         "Overall Acceptance Correction Factors",
                         100, 0.0, 0.18)  # adjust binning & range as needed

    # Keep references alive to avoid garbage collection
    canvas_list = []
    hist_list = []

    for q2 in range(acceptance.shape[0]):
        for xb in range(acceptance.shape[1]):
            for t in range(acceptance.shape[2]):

                # acceptance vs phi
                h_acc = ROOT.TH1D(
                    f"h_acc_q{q2}_xb{xb}_t{t}",
                    f"Acceptance Q2[{Q2_edges[q2]:.2f}-{Q2_edges[q2+1]:.2f}], "
                    f"Xb[{Xb_edges[xb]:.2f}-{Xb_edges[xb+1]:.2f}], "
                    f"-t[{t_edges[t]:.2f}-{t_edges[t+1]:.2f}]",
                    len(phi_edges)-1, phi_edges
                )
                h_acc.SetFillColor(ROOT.kGreen-9)
                h_acc.SetLineColor(ROOT.kGreen+2)
                h_acc.SetFillStyle(1001)

                # generated vs phi
                h_gen = ROOT.TH1D(f"h_gen_q{q2}_xb{xb}_t{t}",
                                  f"Generated Counts;#phi (deg);Counts",
                                  len(phi_edges)-1, phi_edges)
                h_gen.SetLineColor(ROOT.kBlue+2)
                h_gen.SetFillColor(ROOT.kBlue-9)
                h_gen.SetFillStyle(1001)

                # reconstructed vs phi
                h_rec = ROOT.TH1D(f"h_rec_q{q2}_xb{xb}_t{t}",
                                  f"Reconstructed Counts;#phi (deg);Counts",
                                  len(phi_edges)-1, phi_edges)
                h_rec.SetLineColor(ROOT.kRed+2)
                h_rec.SetFillColor(ROOT.kRed-9)
                h_rec.SetFillStyle(1001)

                # fill bins
                for phi in range(len(phi_edges)-1):
                    val = acceptance[q2, xb, t, phi]
                    h_acc.SetBinContent(phi+1, val)
                    if (not np.isnan(val)) and val > 0:
                        h_global.Fill(val)

                    h_gen.SetBinContent(phi+1, gen_counts[q2, xb, t, phi])
                    h_rec.SetBinContent(phi+1, rec_counts[q2, xb, t, phi])

                # only process non-empty histograms
                if h_gen.Integral() > 0 or h_rec.Integral() > 0:

                    # --- overlay canvas with frame ---
                    c = ROOT.TCanvas(f"c_genrec_q{q2}_xb{xb}_t{t}",
                                     "Generated vs Reconstructed", 800, 600)
                    c.SetLogy() 
                    # create a frame histogram to set axes and title
                    frame = ROOT.TH1D(f"frame_q{q2}_xb{xb}_t{t}",
                                      "Generated vs Reconstructed;#phi (deg);Counts",
                                      len(phi_edges)-1, phi_edges)
                    # set y-range slightly above max of gen/rec for visibility
                    ymax = max(h_gen.GetMaximum(), h_rec.GetMaximum())*1.2
                    frame.SetMaximum(ymax)
                    frame.SetMinimum(1e-1)
                    frame.Draw()  # draw axes/frame

                    # draw the histograms on top
                    h_gen.Draw("HIST SAME")
                    h_rec.Draw("HIST SAME")

                    # legend
                    leg = ROOT.TLegend(0.65, 0.75, 0.9, 0.9)
                    leg.AddEntry(h_gen, "Generated", "f")
                    leg.AddEntry(h_rec, "Reconstructed", "f")
                    leg.Draw()

                    c.Update()
                    c.Write()

                    h_acc.Write()

                    # keep references alive
                    canvas_list.append(c)
                    hist_list.extend([h_gen, h_rec, h_acc, frame])

    # global acceptance histogram
    h_global.Write()
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

def luminosity(q_beam, l_target=TARGET_LENGTH, rho_target=RHO_LH2, molar_mass=2.016):
    # returns luminosity in inverse femtobarns, provided input beam charge in C
    return N_A * l_target * rho_target / molar_mass * q_beam / ELECTRON_CHARGE * 1e-39  # fb^-1

### ----------------- ARGS ----------------- ###
parser = argparse.ArgumentParser(description="Compute XSection from input ROOT file.")
parser.add_argument("input_file", type=str)
parser.add_argument("E", type=float, help="Beam energy")
parser.add_argument("--gemc", type=str, help="Optional GEMC file for acceptance corrections")
parser.add_argument("-a","--adaptive", action="store_true")
parser.add_argument("-c6", "--clas6", action="store_true")
parser.add_argument("-c12", "--clas12", action="store_true")
parser.add_argument("-m", "--manual", action="store_true")
parser.add_argument("--local_cuts", action="store_true", help="Perform exclusivity cuts locally per bin instead of global cuts (default)")
parser.add_argument("-v", "--verbose", action="store_true", help="Verbose, diagnostics for all kinematic bins.")
args = parser.parse_args()

# Open the log file once
if args.verbose:
    log_file = open("xsec_log.txt", "w")
    sys.stdout = log_file


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
    Xb_edges = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.38, 0.48, 0.58, 0.7])
    t_edges  = np.array([0.09, 0.15, 0.2, 0.3, 0.4, 0.6, 1.0, 1.5, 2.0])
    phi_edges = np.linspace(0, 360, 21)
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
print(f"Total REC events in designated kinematic range: {np.sum(in_range)}")

bins_map = (len(Q2_edges)-1, len(Xb_edges)-1, len(t_edges)-1, len(phi_edges)-1, 2)

### ------ Sideband subtraction parameters:
n_sigma_signal = 3
n_sigma_sb_min = 3
n_sigma_sb_max = 5

# Flags to optionally toggle cuts
use_global_cuts = not args.local_cuts  # global cuts are default
apply_cuts = True
apply_sideband = True

### -------- Global ex. cuts (per proton topology) --------- ###
global_cuts = {}
global_mgg = {}  # store global m_gg fit per topology
if apply_cuts and use_global_cuts:
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

            # sequentially apply exclusivity cut
            surv_mask[surv_mask] &= np.abs(vals - mu) <= n_sigma_signal*sigma
            if args.verbose:
                print(f"Var={var}, pIdx={p}: mu={mu:.5f}, sigma={sigma:.5f}, survived {np.sum(surv_mask)}")

        # after exclusivity cuts, determine global mgg mean & sigma
        if apply_sideband:
            mgg_vals = df_numpy["eppi0.m_gg"][surv_mask]
            if len(mgg_vals) > 0:
                mu_mgg, sigma_mgg = fit_exclusive(mgg_vals, name="m_gg_global", nbins=200)
                global_mgg[p] = (mu_mgg, sigma_mgg)
                if args.verbose:
                    print(f"Global m_gg, pIdx={p}: mu={mu_mgg:.5f}, sigma={sigma_mgg:.5f}")

### ------ create 4D arrays + boolean mask for final survivors:
side_sub_yield4D = np.zeros(bins_map[:4], dtype=float) 
errs4D = np.zeros(bins_map[:4], dtype=float)
xsec4D = np.zeros(bins_map[:4], dtype=float) 
surv_mask_all = np.zeros(len(df_numpy["dis.Q2"]), dtype=bool)

BEAM_Q = summary_tree.TotalCharge * 1e-9 # C
# print("Beam Q [C] = ", BEAM_Q)
LUM_INT = luminosity(BEAM_Q) # fb^-1
#LUM_INT = 20 # fb^-1
print("Integrated luminosity = ", LUM_INT)
BR = 0.988

### ----------------- Acceptance from GEMC: ----------------- ###
acceptance = None
if args.gemc:
    cols_gemc = ["gen_dis.Q2","gen_dis.Xb","gen_eppi0.t","gen_eppi0.trentoPhi", "dis.Q2","dis.Xb","eppi0.t","eppi0.trentoPhi"]
    df_numpy_gemc = load_gemc(args.gemc, cols_gemc)
    valid_mask = is_valid_rec(df_numpy_gemc)
    print(f"Non-sentinel GEMC reconstructed events: {np.sum(valid_mask)}")
    acceptance, gen_counts, rec_counts = compute_acceptance(df_numpy_gemc, Q2_edges, Xb_edges, t_edges, phi_edges, xsec4D.shape, valid_mask)
    save_acceptance_maps(acceptance, gen_counts, rec_counts, Q2_edges, Xb_edges, t_edges, phi_edges)

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
                    if apply_cuts:
                        for var in exclusive_vars:
                            if args.verbose and np.sum(surv_mask) > 0:
                                print(f"Initiating {var} cuts with {np.sum(surv_mask)} events.")
                            if use_global_cuts:
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
                    if apply_sideband:
                        mgg_vals = df_numpy["eppi0.m_gg"][surv_mask]

                        if use_global_cuts:
                            if p not in global_mgg:
                                continue
                            mu_mgg, sigma_mgg = global_mgg[p]
                        else:
                            if len(vals) == 0:
                                continue
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
                        alpha = w_sig / w_sb if w_sb > 0 else 0.0
                        n_bkg = n_sb * alpha
                        if args.verbose and n_sig > 0:
                            print(f"Bin [Q2: {q2}, Xb: {xb}, t: {t}, phi: {phi}] contains {n_sig} signal events and {n_bkg} background events.")
                        bin_variance += n_sig + (alpha**2) * n_sb
                        bin_yield += max(n_sig - n_bkg, 0.0)

                    else:
                        # No sideband subtraction
                        bin_yield += np.sum(surv_mask)
                        bin_variance += np.sum(surv_mask)

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

core_vars = {
    "Q2":   ("dis.Q2",  "Q^{2} [GeV^{2}]"),
    "xB":   ("dis.Xb",  "x_{B}"),
    "t":    ("eppi0.t", "-t [GeV^{2}]"),
    "phi":  ("eppi0.trentoPhi", "#phi_{Trento} [deg]"),
}

# Fill histograms for core kinematic variables
for key, (branch, title) in core_vars.items():
    vals = df_numpy[branch][surviving_indices]
    if len(vals) == 0:
        continue

    # Special case for phi: map into [0,360)
    if key == "phi":
        vals = (vals + 2*np.pi) % (2*np.pi) * 180.0/np.pi

    h = ROOT.TH1D(
        f"h_{key}_all",
        f"{title};{title};Events",
        nbins_default,
        np.min(vals),
        np.max(vals)
    )
    for v in vals:
        h.Fill(float(v))
    h.Write()
    hist_dict[key] = h

# --- Extra histograms for variables already in ex_vars ---
for var in ex_vars:
    vals = df_numpy[f"eppi0.{var}"][surviving_indices]
    if len(vals) == 0:
        continue
    h = ROOT.TH1D(f"h_{var}_all", f"h_{var}_all", nbins_default, np.min(vals), np.max(vals))
    for v in vals:
        h.Fill(float(v))
    h.Write()
    hist_dict[var] = h

# --- 2D coverage: Q2 vs xB ---
h_Q2_Xb = ROOT.TH2D(
    "h_Q2_Xb", "Q^{2} vs x_{B} Coverage; x_{B}; Q^{2} [GeV^{2}]; Events",
    len(Xb_edges)-1, Xb_edges,
    len(Q2_edges)-1, Q2_edges
)
for q2, xb in zip(df_numpy["dis.Q2"][surviving_indices], df_numpy["dis.Xb"][surviving_indices]):
    h_Q2_Xb.Fill(xb, q2)
h_Q2_Xb.Write()

tree_raw.AutoSave("FlushBaskets")
out_file.Close()
print("Saved exclusive histograms to survival_summary.root")

### ----------------- Project 4D histogram to phi ----------------- ###
out_file = ROOT.TFile("phi_xsec.root", "RECREATE")
Ebeam = args.E

phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
phi_widths = 0.5 * (phi_edges[1:] - phi_edges[:-1])

graph_list = []
fit_list = []

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
            t_center  = np.mean(df_numpy["eppi0.t"][mask_bin])

            # compute virtual photon flux
            gamma = gamma_flux(Q2_center, Xb_center, Ebeam)
            if gamma <= 0:
                continue

            # divide by flux to get reduced cross section
            xsec_red = phi_yields / gamma
            xsec_red_err = phi_errors / gamma

            # skip if xsec_red has zero variance
            if np.allclose(np.std(xsec_red), 0):
                print(f"Skipping bin Q2={q2}, Xb={xb}, t={t} due to zero variance in xsec_red")
                continue

            n_points = len(phi_centers)
            gr_name = f"gr_phi_q{q2}_xb{xb}_t{t}"
            gr_title = (
                f"<Q2>={Q2_center:.2f}, "
                f"<Xb>={Xb_center:.2f}, "
                f"<-t>={t_center:.2f}; "
                f"#phi [deg]; #pi^0 Reduced Cross Section [fb/GeV^2]"
            )

            # --- TGraphAsymmErrors for scatter points with error bars ---
            gr = ROOT.TGraphAsymmErrors(n_points)
            for i in range(n_points):
                gr.SetPoint(i, phi_centers[i], xsec_red[i])
                gr.SetPointError(i, phi_widths[i], phi_widths[i], xsec_red_err[i], xsec_red_err[i])

            # gr.SetMarkerStyle(20)  # solid circle
            # gr.SetMarkerSize(1.0)
            gr.SetName(gr_name)
            gr.SetTitle(gr_title)

            # --- Fit: A + B cos(phi) + C cos(2 phi) ---
            fit_func = ROOT.TF1("fit_phi", "[0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad())", 0, 360)
            fit_func.SetParameters(np.mean(xsec_red), 0.1*np.mean(xsec_red), 0.1*np.mean(xsec_red))
            
            # safer initial parameters
            mean_val = np.mean(xsec_red)
            fit_func.SetParameters(max(mean_val, 1e-6), 0.1*abs(mean_val)+1e-6, 0.1*abs(mean_val)+1e-6)
            fit_func.SetParLimits(0, -2, 2*max(mean_val, 1e-6))
            fit_func.SetParLimits(1, -2*max(mean_val, 1e-6), 2*max(mean_val, 1e-6))
            fit_func.SetParLimits(2, -2*max(mean_val, 1e-6), 2*max(mean_val, 1e-6))

            # fit quietly and robustly
            fit_status = gr.Fit(fit_func, "QR")  # Q=quiet, R=use fit range (?)
            if fit_status != 0:
                print(f"Fit warning in bin Q2={q2}, Xb={xb}, t={t}")

            gr.Draw("AP")         # scatter points with errors
            fit_func.Draw("same") # smooth fit on top
            gr.GetXaxis().SetRangeUser(0, 360)

            q2_mean_param = ROOT.TParameter("double")("q2_mean", Q2_center)
            xb_mean_param = ROOT.TParameter("double")("xb_mean", Xb_center)
            t_mean_param  = ROOT.TParameter("double")("t_mean", t_center)
            gr.GetListOfFunctions().Add(q2_mean_param)
            gr.GetListOfFunctions().Add(xb_mean_param)
            gr.GetListOfFunctions().Add(t_mean_param)
            
            gr.Write()
            graph_list.append(gr)      
            fit_list.append(fit_func)  

out_file.Close()

if args.verbose:
    sys.stdout = sys.__stdout__   # restore original stdout
    log_file.close()
    print("Verbose output saved to xsec_log.txt")

print("Saved phi-projected reduced cross sections [nb/GeV^2] to phi_xsec.root")

