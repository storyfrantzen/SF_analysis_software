# type: ignore
import sys
import argparse
import ROOT
import numpy as np
import math

ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

# constants:
ELECTRON_CHARGE = 1.602e-19 # Coulombs
TARGET_LENGTH = 5.00 # cm
RHO_LH2 = 0.071 # g / cm3
N_A = 6.02214e23
ALPHA_EM = 1/137.035999084
PROTON_MASS = 0.9382720813  # proton mass [GeV]
BR = 0.988 # branching ratio for pi0 --> 2g
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
def get_adaptive_edges(tree, var_name, n_bins, min_val, max_val, fine_bins=1000):
    hist_name = f"h_{var_name}"
    hist = ROOT.TH1D(hist_name, hist_name, fine_bins, min_val, max_val)
    tree.Draw(f"{var_name} >> {hist_name}", "", "goff")
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

def load_np_from_file(file, cols):
    """Load ROOT file into numpy dict."""
    f = ROOT.TFile(file, "READ")
    tree = f.Get("Events")
    df = ROOT.RDataFrame(tree)
    df_np = df.AsNumpy(columns=cols)
    for k in cols:
        df_np[k] = coerce_scalar_column(df_np[k])
    return df_np

def no_sentinels_mask(df_np):
    mask = np.ones(len(df_np["dis.Q2"]), dtype=bool)
    for k, arr in df_np.items():
        if k.startswith("eppi0.") or k.startswith("dis."):
            if arr.dtype.kind in ("i", "u"):  # integers
                mask &= (arr != -9999)
            else:  # floats
                mask &= ~np.isnan(arr)
    return mask

def fit_exclusive(vals, var_name=None, fit_func=None, fit_range=None, p0=None, nbins=None, fallback=False, min_entries=10):
    """
    Fit a 1D distribution and return the mean and sigma.
    
    Parameters
    ----------
    vals : array-like
        Data points to fit.
    var_name : str
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
    vals = np.asarray(vals, float)
    vals = vals[~np.isnan(vals)]
    n = len(vals)
    if n == 0:
        print("Fit error: no data found!")
        return np.nan, np.nan

    if nbins is None:
        nbins = max(30, min(200, int(max(10, n // 2))))
        print(f"Defaulting to {nbins} for {var_name}.")
    if n < min_entries:
        print(f"Fit error: fewer than {min_entries} events, cannot perform fit!")
        if fallback:
            mu, sigma = np.mean(vals), np.std(vals, ddof=1) if n > 1 else np.nan
            return mu, sigma
        return np.nan, np.nan

    # Configure histograms for fitting:

    min_data, max_data = float(np.min(vals)), float(np.max(vals))

    hist_name = f"hfit_{np.random.randint(1e9)}"
    hist = ROOT.TH1D(hist_name, hist_name, nbins, min_data, max_data)
    for v in vals: 
        hist.Fill(v)

    # Determine default fit settings from dict
    default_info = DEFAULT_FITS.get(var_name, {})

    # Set fit range
    fit_range_tuple  = fit_range or default_info.get("range", (min_data, max_data))

    if fit_range_tuple[0] < min_data:
        print("Fit range warning: setting min_fit to min_data.")
    if fit_range_tuple[1] > max_data:
        print("Fit range warning: setting max_fit to max_data.")
    min_fit = max(fit_range_tuple[0], min_data)
    max_fit = min(fit_range_tuple[1], max_data)
    if min_fit >= max_fit:
        print("Fit range error: fit range collapsed!")
        if fallback:
            return np.mean(vals), np.std(vals, ddof=1)
        return np.nan, np.nan
    fit_range_tuple = (min_fit, max_fit)

    # Set fit function string / TF1 object
    if fit_func is None:
        func_str = default_info.get("fit_func", "gaus(0)+pol1(3)")
        tf1_obj = None
    else:
        if isinstance(fit_func, ROOT.TF1):
            func_str = fit_func.GetTitle()
            tf1_obj = fit_func
        else:  # treat as formula string
            func_str = fit_func
            tf1_obj = None

    if tf1_obj:
        f = tf1_obj.Clone(f"f_{hist_name}")
        f.SetRange(min_fit, max_fit)
    else:
        f = ROOT.TF1(f"f_{hist_name}", func_str, min_fit, max_fit)

    # Set initial parameters
    p0_tuple = p0 or default_info.get("p0", None)

    if p0_tuple is None:
        amp   = hist.GetMaximum()
        peak  = hist.GetBinCenter(hist.GetMaximumBin())
        width = hist.GetStdDev() if hist.GetStdDev() > 0 else 0.05
        if "crystalball" in func_str.lower():
            # amp, μ, σ, alpha (-1 * where the tail begins), n (tail power)
            p0_tuple = (amp, peak, width, -1.5, 4)
        else: 
            # gaussian + poly2/1 etc.
            # (amp, mu, sigma, <poly params...>)
            if "pol1" in func_str:
                p0_tuple = (amp, peak, width, 0.0, 0.0)
            elif "pol2" in func_str:
                p0_tuple = (amp, peak, width, 0.0, 0.0, 0.0)
            else:
                p0_tuple = (amp, peak, width)

    for i, val in enumerate(p0_tuple):
        f.SetParameter(i, val)

    # Perform fit and extract results
    status = hist.Fit(f, "RQ0")
    fit_failed = (status != 0)

    if fit_failed:
        if fallback:
            mu = np.mean(vals)
            sigma = np.std(vals, ddof=1)
        else:
            "Fit error: fit failed! Returning np.nan for mu, sigma."
            mu = sigma = np.nan
    else:
        try:
            mu = f.GetParameter(1)
            sigma = abs(f.GetParameter(2))
        except Exception:
            if fallback:
                mu = np.mean(vals)
                sigma = np.std(vals, ddof=1)
            else:
                mu = sigma = np.nan

    # cleanup
    del hist
    return mu, sigma

def compute_global_cuts(df_np, in_range_mask, pIdx, exclusive_vars, n_sigma_signal=3, apply_sideband=False, verbose=False, label="DATA"):
    """
    Compute global exclusivity cuts (and optional m_gg sideband cuts)
    for either DATA or GEMC.

    Returns:
        global_cut_params: dict mapping (var, pIdx) -> (mu, sigma)
        global_exclusivity_mask: boolean mask of events surviving ALL cuts
    """

    N = len(df_np["dis.Q2"])
    global_exclusivity_mask = np.zeros(N, dtype=bool)
    global_cut_params = {}

    print(f">>> Deriving global exclusivity cuts for {label}...")

    for p in range(2):  # pIdx = 0 (FD), 1 (CD)
        mask = (in_range_mask) & (pIdx == p)

        for var in exclusive_vars:
            vals = df_np[f"eppi0.{var}"][mask]

            if len(vals) == 0:
                continue

            # Fit exclusivity variable globally
            mu, sigma = fit_exclusive(vals, var_name=var)
            global_cut_params[(var, p)] = (mu, sigma)

            # Apply the cut
            mask[mask] &= np.abs(vals - mu) <= n_sigma_signal * sigma
            global_exclusivity_mask |= mask

            if verbose:
                print(f"{label}: Global {var}, pIdx={p}: "
                      f"mu={mu:.5f}, sigma={sigma:.5f}, survived={np.sum(mask)}")

        # Optional m_gg global sideband fits
        if apply_sideband:
            vals = df_np["eppi0.m_gg"][mask]
            if len(vals) > 0:
                mu, sigma = fit_exclusive(vals, var_name="m_gg", nbins=200)
                global_cut_params[("m_gg", p)] = (mu, sigma)
                if verbose:
                    print(f"{label}: Global m_gg, pIdx={p}: mu={mu:.5f}, sigma={sigma:.5f}")

    return global_cut_params, global_exclusivity_mask

def compute_acceptance(df_np_gemc, bin_edges, shape, valid_mask=None, min_acc=0.005, cfg={}):
    """Return acceptance[Q2, Xb, t, phi] from GEMC numpy arrays."""

    # --- unpack required bin edges ---
    Q2_edges, Xb_edges, t_edges, phi_edges = bin_edges

    # --- unpack optional config ---
    use_global_cuts = cfg.get("use_global_cuts", False)
    apply_sideband  = cfg.get("apply_sideband", False)
    verbose         = cfg.get("verbose", False)
    n_sigma_signal  = cfg.get("n_sigma_signal", 3)
    exclusive_vars  = cfg.get("exclusive_vars", [])

    # output arrays
    gen_counts = np.zeros(shape, dtype=float)
    rec_counts = np.zeros(shape, dtype=float)

    # indices for generated
    iQ2_gen  = get_bin_indices(df_np_gemc["gen_dis.Q2"], Q2_edges)
    iXb_gen  = get_bin_indices(df_np_gemc["gen_dis.Xb"], Xb_edges)
    it_gen   = get_bin_indices(df_np_gemc["gen_eppi0.t"],  t_edges)
    iphi_gen = get_bin_indices((df_np_gemc["gen_eppi0.trentoPhi"]+2*np.pi)%(2*np.pi)*180/np.pi, phi_edges)
    valid_gen = (iQ2_gen>=0)&(iXb_gen>=0)&(it_gen>=0)&(iphi_gen>=0)
    print(f"Total GEN gemc events in kinematic binning range: {np.sum(valid_gen)}")
    
    # indices for reconstructed
    iQ2_rec  = get_bin_indices(df_np_gemc["dis.Q2"], Q2_edges)
    iXb_rec  = get_bin_indices(df_np_gemc["dis.Xb"], Xb_edges)
    it_rec   = get_bin_indices(df_np_gemc["eppi0.t"],  t_edges)
    iphi_rec = get_bin_indices((df_np_gemc["eppi0.trentoPhi"]+2*np.pi)%(2*np.pi)*180/np.pi, phi_edges)
    valid_rec = (iQ2_rec>=0) & (iXb_rec>=0) & (it_rec>=0) & (iphi_rec>=0)
    print(f"Total REC gemc events in kinematic binning range: {np.sum(valid_rec)}")
    if valid_mask is not None:
        valid_rec &= valid_mask
    pIdx_rec = np.where(np.isclose(df_np_gemc["p.det"], 1.0), 0, 1) # logic might need to be checked...

    global_mask_gemc = np.ones_like(valid_rec, dtype=bool)
    
    if use_global_cuts:
        global_cut_params_gemc, global_mask_gemc = compute_global_cuts(
            df_np=df_np_gemc,
            in_range_mask=valid_rec,
            pIdx=pIdx_rec,
            exclusive_vars=exclusive_vars,
            apply_sideband=apply_sideband,
            n_sigma_signal=n_sigma_signal,
            verbose=verbose,
            label="GEMC"
        )
    
    valid_rec &= global_mask_gemc

    np.add.at(gen_counts, (iQ2_gen[valid_gen], iXb_gen[valid_gen], it_gen[valid_gen], iphi_gen[valid_gen]), 1)
    np.add.at(rec_counts, (iQ2_rec[valid_rec], iXb_rec[valid_rec], it_rec[valid_rec], iphi_rec[valid_rec]), 1)

    acceptance = np.divide(rec_counts, gen_counts, out=np.zeros_like(rec_counts), where=gen_counts>0)

    # --- mask bins below threshold ---
    acceptance_masked = np.where(acceptance >= min_acc, acceptance, np.nan)

    return acceptance_masked, gen_counts, rec_counts

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
                    if not np.isnan(val):
                        h_acc.SetBinContent(phi+1, val)
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
                    leg.AddEntry(h_gen, f"Generated (N={int(h_gen.Integral())})", "f")
                    leg.AddEntry(h_rec, f"Reconstructed (N={int(h_rec.Integral())})", "f")
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

### ------ Sideband subtraction parameters:
n_sigma_signal = 3
n_sigma_sb_min = 3
n_sigma_sb_max = 5

# Flags to toggle cuts
use_global_cuts = not args.local_cuts  
apply_cuts = True
apply_sideband = False

# Exclusivity variables:
ex_vars = ["m2_epX", "m_eggX", "E_miss", "pi0_thetaX", "m_gg"]
#exclusive_vars = ["m2_epX", "m_eggX", "E_miss", "pi0_thetaX"] # not including m_gg
exclusive_vars = ["m2_epX", "m_eggX", "E_miss", "pi0_thetaX", "m_gg"] 

### --------------- Initialize bins --------------- ###
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
    raise RuntimeError("Must specify binning scheme, e.g., --adaptive, --clas12, --clas6, or --manual")

bins_map = (len(Q2_edges)-1, len(Xb_edges)-1, len(t_edges)-1, len(phi_edges)-1, 2)

### ------ create 4D arrays + boolean mask for final survivors:
side_sub_yield4D = np.zeros(bins_map[:4], dtype=float) 
xsec4D = np.zeros(bins_map[:4], dtype=float) 
errs4D = np.zeros(bins_map[:4], dtype=float)

### --------------- Read eppi0REC ROOT file --------------- ###
f = ROOT.TFile(args.input_file, "READ")
events_tree = f.Get("Events")
summary_tree = f.Get("Summary")
summary_tree.GetEntry(0)
BEAM_Q = summary_tree.TotalCharge * 1e-9 # C
LUM_INT = luminosity(BEAM_Q) # fb^-1
print("Beam Q [C] = ", BEAM_Q, " ; Integrated Luminosity [fb^-1] = ", LUM_INT)

cols_rec = ["dis.Q2","dis.Xb","eppi0.t","eppi0.trentoPhi","p.det","eppi0.pi0_thetaX","eppi0.m2_epX","eppi0.m_eggX","eppi0.E_miss","eppi0.m_gg"]
df_np_data = load_np_from_file(args.input_file, cols_rec)

iQ2_data   = get_bin_indices(df_np_data["dis.Q2"], Q2_edges)
iXb_data   = get_bin_indices(df_np_data["dis.Xb"], Xb_edges)
it_data    = get_bin_indices(df_np_data["eppi0.t"], t_edges)
iphi_data  = get_bin_indices((df_np_data["eppi0.trentoPhi"] + 2*np.pi) % (2*np.pi) * 180.0/np.pi, phi_edges)
pIdx_data  = np.where(np.isclose(df_np_data["p.det"], 1.0), 0, 1).astype(int)

in_range_data = (iQ2_data>=0) & (iXb_data>=0) & (it_data>=0) & (iphi_data>=0)
print(f"Total REC data events in kinematic binning range: {np.sum(in_range_data)}")

### ----------------- Acceptance from GEMC: ----------------- ###
acceptance = None
if args.gemc:
    cols_gemc = ["gen_dis.Q2","gen_dis.Xb","gen_eppi0.t","gen_eppi0.trentoPhi"] + list(cols_rec)
    df_np_gemc = load_np_from_file(args.gemc, cols_gemc)
    total_gen_events = len(df_np_gemc["gen_dis.Q2"])
    print(f"Total generated GEMC events: {total_gen_events}")
    valid_mask = no_sentinels_mask(df_np_gemc)
    cfg = {}
    if apply_cuts:
        cfg = {"use_global_cuts": use_global_cuts, 
               "exclusive_vars": exclusive_vars, 
               "apply_sideband": apply_sideband, 
               "n_sigma_signal": n_sigma_signal, 
               "verbose": args.verbose}
    acceptance, gen_counts, rec_counts = compute_acceptance(df_np_gemc, (Q2_edges, Xb_edges, t_edges, phi_edges), xsec4D.shape, valid_mask, cfg=cfg)
    save_acceptance_maps(acceptance, gen_counts, rec_counts, Q2_edges, Xb_edges, t_edges, phi_edges)


### -------- Global ex. cuts (per proton topology) --------- ###
if apply_cuts and use_global_cuts:
    global_cut_params_data, global_exclusivity_mask_data = compute_global_cuts(
        df_np=df_np_data,
        in_range_mask=in_range_data,
        pIdx=pIdx_data,
        exclusive_vars=exclusive_vars,
        apply_sideband=apply_sideband,
        n_sigma_signal=n_sigma_signal,   # can rework these kwargs with cfg from above
        verbose=args.verbose,
        label="DATA"
    )

### ----------- Create tree to track raw counts ----------- ###
out_file = ROOT.TFile("survival_summary.root", "RECREATE")
tree_raw = ROOT.TTree("raw_counts", "Raw event counts per bin")
tree_raw.SetAutoFlush(10000)  # flush every 10k entries

q2b   = np.zeros(1, dtype=int)
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

### ----------- Helpers for bins near acceptance boundary ----------- ###
Ebeam = args.E
W_min = 2.0   # GeV
Q2_min = 1.0  # GeV^2

# numerical integration settings
N_rie = 100             # number of samples per xB bin (Riemann sum)
eprime_scale = 0.81     # 1 - min(p_e')/E_beam, ONLY A PLACEHOLDER BASED ON RGK VALUES CURRENTLY
Q2_left_slope = 2 * PROTON_MASS * Ebeam * eprime_scale

### ---------------------------------- XSECTION Workflow ---------------------------------- ###
surv_mask_all = np.zeros(len(df_np_data["dis.Q2"]), dtype=bool)

for q2 in range(bins_map[0]):
    for xb in range(bins_map[1]):
        for t in range(bins_map[2]):
            dt = t_edges[t+1] - t_edges[t]
            for phi in range(bins_map[3]):
                dphi = phi_edges[phi+1] - phi_edges[phi]
        
                mask_bin = (in_range_data) & (iQ2_data==q2) & (iXb_data==xb) & (it_data==t) & (iphi_data==phi)
                num_events = np.sum(mask_bin)
                if num_events > 0:
                    q2b[0], xbb[0], tb[0], phib[0], count[0] = q2, xb, t, phi, num_events
                    tree_raw.Fill()
                
                x_low, x_high = Xb_edges[xb], Xb_edges[xb+1]
                dx = (x_high - x_low) / float(N_rie)
                # midpoints: avoid endpoints to be stable at singularities
                x_midpoints = np.linspace(x_low + dx/2.0, x_high - dx/2.0, N_rie)
                # evaluate physical Q2 bounds at the sampled x points
                q2_phys_low = np.maximum(Q2_min, ((W_min**2 - PROTON_MASS**2) / (1.0/x_midpoints - 1.0)))   # physical lower bound at each x
                q2_phys_high = Q2_left_slope * x_midpoints      # physical upper bound at each x
                # compute overlap of the bin's [q2_low, q2_high] with physical [q2_phys_low, q2_phys_high]
                q2_bin_low, q2_bin_high = Q2_edges[q2], Q2_edges[q2+1]
                # elementwise overlap height: max(0, min(bin_high, phys_high) - max(bin_low, phys_low))
                overlap_low = np.maximum(q2_bin_low, q2_phys_low)
                overlap_high = np.minimum(q2_bin_high, q2_phys_high)
                overlap_height = np.maximum(0.0, overlap_high - overlap_low)  # array of length N_rie
                # approximate area in Q2-xB for this (q2,xb) bin by Riemann sum
                area_Q2_xB = np.sum(overlap_height) * dx   # [GeV^2 * xB]

                # RECTANGULAR BIN DEFAULT:
                # dq2 = Q2_edges[q2+1] - Q2_edges[q2]
                # dxb = Xb_edges[xb+1] - Xb_edges[xb]

                bin_volume = area_Q2_xB * dt * dphi
                if bin_volume == 0: # or not np.any(mask_bin):
                    continue

                bin_yield = 0
                bin_variance = 0

                # apply global cuts once
                if apply_cuts and use_global_cuts:
                    mask_bin &= global_exclusivity_mask_data
                
                # check if we need the proton loop
                use_proton_loop = (apply_cuts and not use_global_cuts) or apply_sideband

                if use_proton_loop:
                    for p in range(2):
                        mask_bin_p = mask_bin & (pIdx_data == p)
                        if not np.any(mask_bin_p):
                            continue
                        if args.verbose:
                            print(f"Bin [Q2:{q2}, Xb:{xb}, t:{t}, phi:{phi}, pIdx:{p}] += {np.sum(mask_bin_p)} events.")

                        local_mask = mask_bin_p.copy()

                        # local cuts
                        if apply_cuts and not use_global_cuts:
                            for var in exclusive_vars:
                                vals = df_np_data[f"eppi0.{var}"][local_mask]
                                mu, sigma = fit_exclusive(vals, name=var)
                                local_mask[local_mask] &= np.abs(vals - mu) <= n_sigma_signal * sigma

                        # sideband subtraction / yield
                        if apply_sideband:
                            mgg_vals = df_np_data["eppi0.m_gg"][local_mask]
                            if use_global_cuts:
                                mu_mgg, sigma_mgg = global_cut_params_data.get(("m_gg", p), (0, 0)) # Note default...
                            else:
                                mu_mgg, sigma_mgg = fit_exclusive(mgg_vals, name="m_gg")

                            sig_mask = np.abs(mgg_vals - mu_mgg) < n_sigma_signal*sigma_mgg
                            sb_mask  = ((mu_mgg - n_sigma_sb_max*sigma_mgg <= mgg_vals) & 
                                        (mgg_vals < mu_mgg - n_sigma_sb_min*sigma_mgg)) | \
                                    ((mu_mgg + n_sigma_sb_min*sigma_mgg < mgg_vals) & 
                                        (mgg_vals <= mu_mgg + n_sigma_sb_max*sigma_mgg))

                            n_sig = np.sum(sig_mask)
                            n_sb  = np.sum(sb_mask)
                            alpha = (2*n_sigma_signal*sigma_mgg) / (2*(n_sigma_sb_max-n_sigma_sb_min)*sigma_mgg) if n_sb > 0 else 0
                            n_bkg = n_sb * alpha
                            bin_variance += n_sig + (alpha**2) * n_sb
                            bin_yield += n_sig - n_bkg
                        else:
                            bin_yield += np.sum(local_mask)
                            bin_variance += np.sum(local_mask) # Need to consider yield == 0 more carefully! 

                        # mark surviving events
                        surv_mask_all[local_mask] = True

                else:
                    bin_yield = np.sum(mask_bin)
                    bin_variance = np.sum(mask_bin)
                    surv_mask_all[mask_bin] = True

                # fill 4D arrays
                if bin_yield == 0 and bin_variance == 0: # need to consider this part more carefully
                    side_sub_yield4D[q2, xb, t, phi] = np.nan
                    errs4D[q2, xb, t, phi] = np.nan
                    xsec4D[q2, xb, t, phi] = np.nan
                else:
                    side_sub_yield4D[q2, xb, t, phi] = bin_yield
                    errs4D[q2, xb, t, phi] = np.sqrt(bin_variance) / bin_volume / LUM_INT / BR
                    xsec4D[q2, xb, t, phi] = bin_yield / bin_volume / LUM_INT / BR

if acceptance is not None:
    xsec4D /= acceptance
    errs4D /= acceptance
    acc_err = np.sqrt(acceptance * (1 - acceptance) / gen_counts)
    rel_err_data = errs4D / xsec4D  # δN / N
    rel_err_acc  = acc_err / acceptance 
    rel_err_total = np.sqrt(rel_err_data**2 + rel_err_acc**2)
    errs4D = xsec4D * rel_err_total

surviving_indices = np.flatnonzero(surv_mask_all)
print(f"Total surviving events: {len(surviving_indices)}")
if apply_sideband:
    print(f"Total yield after 3σ cuts + sideband subtraction: {np.nansum(side_sub_yield4D):.3f}")

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
    vals = df_np_data[branch][surviving_indices]
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
    vals = df_np_data[f"eppi0.{var}"][surviving_indices]
    if len(vals) == 0:
        continue
    h = ROOT.TH1D(f"h_{var}_all", f"h_{var}_all", nbins_default, np.min(vals), np.max(vals))
    for v in vals:
        h.Fill(float(v))
    h.Write()
    hist_dict[var] = h

# --- 2D coverage: Q2 vs xB ---
# xB points for drawing
xB_vals = np.linspace(min(Xb_edges), max(Xb_edges), 200)

# Compute xB interval between left and right kinematic boundaries
xB_min_lower = Q2_min / Q2_left_slope
xB_max_lower = Q2_min / ( (W_min**2 - PROTON_MASS**2) + Q2_min )
xB_lower = xB_vals[(xB_vals >= xB_min_lower) & (xB_vals <= xB_max_lower)]
q2_lower_trimmed = np.full_like(xB_lower, Q2_min)
xB_higher = xB_vals[xB_vals >= xB_max_lower]

# --- Left boundary line ---
q2_left = Q2_left_slope * xB_vals[xB_vals >= xB_min_lower]
gr_left = ROOT.TGraph(len(xB_vals[xB_vals >= xB_min_lower]), xB_vals[xB_vals >= xB_min_lower], q2_left)
gr_left.SetLineColor(ROOT.kMagenta)
gr_left.SetLineWidth(2)

# --- Lower horizontal line ---
gr_lower = ROOT.TGraph(len(xB_lower), xB_lower, q2_lower_trimmed)
gr_lower.SetLineColor(ROOT.kRed)
gr_lower.SetLineWidth(2)

# --- Right boundary curve ---
q2_right = (W_min**2 - PROTON_MASS**2) / (1/xB_higher - 1)
gr_right = ROOT.TGraph(len(xB_higher), xB_higher, q2_right)
gr_right.SetLineColor(ROOT.kCyan)
gr_right.SetLineWidth(2)

h_Q2_Xb_binned = ROOT.TH2D(
    "h_Q2_Xb_binned", "Q^{2} vs x_{B} Coverage; x_{B}; Q^{2} [GeV^{2}]; Events",
    len(Xb_edges)-1, Xb_edges,
    len(Q2_edges)-1, Q2_edges
)
h_Q2_Xb = ROOT.TH2D(
    "h_Q2_Xb", "Q^{2} vs x_{B} Coverage; x_{B}; Q^{2} [GeV^{2}]; Events",
    100, 0.05, 0.8,
    100, 0.0, np.max(df_np_data["dis.Q2"][surviving_indices])
)
for q2, xb in zip(df_np_data["dis.Q2"][surviving_indices], df_np_data["dis.Xb"][surviving_indices]):
    h_Q2_Xb.Fill(xb, q2)
    h_Q2_Xb_binned.Fill(xb, q2)

h_Q2_Xb_binned.SetFillColorAlpha(ROOT.kOrange, 0.4)  # 0.4 = 40% opaque
h_Q2_Xb_binned.SetLineColor(ROOT.kOrange)  # optional outline
h_Q2_Xb_binned.SetMarkerSize(0)  # no markers

c = ROOT.TCanvas("Q2_vs_Xb", "Q2 vs xB Coverage", 400, 600)
h_Q2_Xb.Draw("COLZ")   # Draw 2D histogram with color map
h_Q2_Xb_binned.Draw("BOX SAME")
gr_lower.Draw("L SAME")
gr_left.Draw("L SAME")
gr_right.Draw("L SAME")

c.Write()

tree_raw.AutoSave("FlushBaskets")
out_file.Close()
print("Saved exclusive histograms to survival_summary.root")

### ----------------- Project 4D histogram to phi ----------------- ###
out_file = ROOT.TFile("phi_xsec.root", "RECREATE")

phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
phi_widths = 0.5 * (phi_edges[1:] - phi_edges[:-1])

graph_list = []
fit_list = []

for q2 in range(bins_map[0]):
    for xb in range(bins_map[1]):
        for t in range(bins_map[2]):
            phi_yields = xsec4D[q2, xb, t, :]
            phi_errors = errs4D[q2, xb, t, :]

            if np.nansum(phi_yields) == 0:
                continue

            # surviving event mask for this 4D bin
            mask_bin = (iQ2_data == q2) & (iXb_data == xb) & (it_data == t) & surv_mask_all
            if not np.any(mask_bin):
                continue

            # compute bin centers for Q2 and Xb (arithmetic mean of surviving events)
            Q2_center = np.mean(df_np_data["dis.Q2"][mask_bin])
            Xb_center = np.mean(df_np_data["dis.Xb"][mask_bin])
            t_center  = np.mean(df_np_data["eppi0.t"][mask_bin])

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

            valid_mask = ~np.isnan(xsec_red)

            gr_name = f"gr_phi_q{q2}_xb{xb}_t{t}"
            gr_title = (
                f"<Q2>={Q2_center:.2f}, "
                f"<Xb>={Xb_center:.2f}, "
                f"<-t>={t_center:.2f}; "
                f"#phi [deg]; #pi^0 Reduced Cross Section [fb/GeV^2]"
            )

            # --- TGraphAsymmErrors for scatter points with error bars ---
            gr = ROOT.TGraphAsymmErrors()
            for i, idx in enumerate(np.where(valid_mask)[0]):
                gr.SetPoint(i, phi_centers[idx], xsec_red[idx])
                gr.SetPointError(i, phi_widths[idx], phi_widths[idx], xsec_red_err[idx], xsec_red_err[idx])

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

print("Saved phi-projected reduced cross sections [fb/GeV^2] to phi_xsec.root")