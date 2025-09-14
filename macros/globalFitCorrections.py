# globalFitCorrections.py
# type: ignore
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import argparse
import json

# Load compiled dictionary for branch variables
ROOT.gSystem.Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so")

def build_design_matrix(theta, p, isCD=False, res=0):
    """
    Construct basis functions for a global fit depending on residual type.

    Parameters
    ----------
    theta : array-like
        theta in degrees (or rad if you pass rad consistently; your script uses degrees).
    p : array-like
        momentum in GeV.
    isCD : bool
        If True, return a design matrix tailored for the Central Detector.
    res : int
        Which residual to build for:
            0 -> delta_p
            1 -> delta_theta (degrees)
            2 -> delta_phi   (degrees)

    Returns
    -------
    X : ndarray, shape (n_events, n_features)
        Design matrix (columns = basis functions).
    """

    theta = np.asarray(theta)
    p = np.asarray(p)

    if isCD:
        # central detector: momentum range and behavior differ → pick different basis
        if res == 0:   # delta_p
            cols = [
                np.ones_like(p),
                p**3,                               
                theta,            
                p**2 * np.log(p) 
            ]
        elif res == 1: # delta_theta (deg)
            cols = [
                np.ones_like(p),
                1.0 / p,        
                np.log(p),       
                theta,
                p**3,                      
            ]
        elif res == 2: # delta_phi (deg)
            cols = [
                np.ones_like(p),
                theta / p,
                p**3,
                1.0 / (p**2)
            ]
        else:
            raise ValueError("res must be 0 (dp), 1 (dtheta) or 2 (dphi)")
    else:
        # Forward detector (pFD) or general-purpose basis
        if res == 0:   # delta_p
            cols = [
                np.ones_like(p),   
                np.log(p),   
                theta/p,
                theta**2/p,
                theta/p**2             
            ]
        elif res == 1: # delta_theta (deg)
            cols = [
                np.ones_like(p),
                np.log(p),
                theta,
                theta/p,
                (theta * np.log(theta))/p,
                theta/p**2
            ]
        elif res == 2: # delta_phi (deg)
            cols = [
                np.ones_like(p),
                np.log(p),
                1/p,
                theta,
                theta/p,
                theta/p**2
            ]
        else:
            raise ValueError("res must be 0 (dp), 1 (dtheta) or 2 (dphi)")

    # stack columns into design matrix
    X = np.column_stack(cols)
    return X

def fixed_global_fit(df, isCD=False):
    """Fit Δp, Δθ, Δφ and return coefficient + chi2 dictionary."""
    arrs = df.AsNumpy(["theta_deg", "rec.p", "delta_p", "delta_theta", "delta_phi"])
    theta = arrs["theta_deg"]
    p = arrs["rec.p"]
    dps = {key: arrs[key] for key in ["delta_p", "delta_theta", "delta_phi"]}

    results = {}
    # map residual name -> index used by build_design_matrix
    res_map = {"delta_p": 0, "delta_theta": 1, "delta_phi": 2}

    for key, y in dps.items():
        # build design matrix tailored to this residual
        res_idx = res_map[key]
        X = build_design_matrix(theta, p, isCD=isCD, res=res_idx)

        # solve least squares
        beta, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
        y_pred = X @ beta
        residuals_arr = y - y_pred

        # compute unweighted chi2
        chi2 = np.sum(residuals_arr**2)
        ndof = len(y) - len(beta)
        chi2_ndof = chi2 / ndof if ndof > 0 else np.nan

        print(f"\nFitted coefficients for {key}:")
        for i, c in enumerate(beta):
            print(f"  c[{i}] = {c:.6e}")
        print(f"  ~χ²/dof = {chi2:.3e} / {ndof} = {chi2_ndof:.9f}")

        results[key] = {
            "beta": beta,
            "chi2": chi2,
            "ndof": ndof,
            "chi2_ndof": chi2_ndof,
            "rank": rank,
            "singular_vals": s.tolist()
        }

    return results

def add_profile(ax, x, y, bins, color="w", label="Profile mean", xrange=None):
    """Overlay profile (mean residual in x-bins) with error bars."""
    if xrange is None:
        xrange = (x.min(), x.max())
    counts, bin_edges = np.histogram(x, bins=bins, range=xrange)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    means, errors = [], []

    for i in range(len(bin_edges) - 1):
        mask = (x >= bin_edges[i]) & (x < bin_edges[i+1])
        if np.any(mask):
            vals = y[mask]
            means.append(np.mean(vals))
            errors.append(np.std(vals) / np.sqrt(len(vals)))  # error on mean
        else:
            means.append(np.nan)
            errors.append(np.nan)

    ax.errorbar(bin_centers, means, yerr=errors, fmt='.-', color=color, 
                label=label, markersize=4, linewidth=1)
    ax.legend(fontsize=8, loc="upper right")

def export_coeffs(coeffs, filename="correction_coeffs.json"):
    out = {k: v["beta"].tolist() for k, v in coeffs.items()}
    with open(filename, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Exported coefficients to {filename}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--tree", default="Events")
    parser.add_argument("--max_events", type=int, default=None)
    parser.add_argument("--isCD", action="store_true")
    parser.add_argument("--draw", action="store_true")
    args = parser.parse_args()

    df = ROOT.RDataFrame(args.tree, args.input_file)

    radtoDeg = 180/np.pi
    df_pro = df.Filter("rec.pid == 2212 && rec.charge == 1")
    df_pro = df_pro.Define("theta_deg", f"rec.theta*{radtoDeg}") \
                 .Define("delta_p", "gen.p - rec.p") \
                 .Define("delta_theta", f"(gen.theta - rec.theta)*{radtoDeg}") \
                 .Define("delta_phi", f"TMath::ATan2(TMath::Sin(gen.phi - rec.phi), TMath::Cos(gen.phi - rec.phi))*{radtoDeg}")
    
    if args.max_events:
        df_pro = df_pro.Range(args.max_events)

    # Detector-specific cuts
    theta_range = [18, 40] if not args.isCD else [39, 60]
    p_range = [0.55, 4] if not args.isCD else [0.3, 2.4]
    df_pro = df_pro.Filter(f"theta_deg > {theta_range[0]} && theta_deg < {theta_range[1]}")
    df_pro = df_pro.Filter(f"rec.det == {1 if not args.isCD else 2}")

    # Fit
    coeffs = fixed_global_fit(df_pro, isCD=args.isCD)

    if args.isCD:
        export_coeffs(coeffs, "pCD_correction_coeffs.json")
    else:
        export_coeffs(coeffs, "pFD_correction_coeffs.json")

    if not args.draw:
        return

    # Build theta/p grid for surfaces
    theta_grid = np.linspace(theta_range[0], theta_range[1], 50)
    p_grid = np.linspace(p_range[0], p_range[1], 50)
    THETA, P = np.meshgrid(theta_grid, p_grid)
    Xgrid_p = build_design_matrix(THETA.ravel(), P.ravel(), isCD=args.isCD, res=0)
    Xgrid_theta = build_design_matrix(THETA.ravel(), P.ravel(), isCD=args.isCD, res=1)
    Xgrid_phi = build_design_matrix(THETA.ravel(), P.ravel(), isCD=args.isCD, res=2)

    delta_p_grid     = (Xgrid_p @ coeffs['delta_p']['beta']).reshape(THETA.shape)
    delta_theta_grid = (Xgrid_theta @ coeffs['delta_theta']['beta']).reshape(THETA.shape)
    delta_phi_grid   = (Xgrid_phi @ coeffs['delta_phi']['beta']).reshape(THETA.shape)

    # Event arrays
    arrs = df_pro.AsNumpy(["theta_deg","rec.p","delta_p","delta_theta","delta_phi"])
    theta = arrs["theta_deg"]
    p = arrs["rec.p"]
    delta_p_orig     = arrs["delta_p"]
    delta_theta_orig = arrs["delta_theta"]
    delta_phi_orig   = arrs["delta_phi"]

    # Corrections & residuals
    X_p = build_design_matrix(theta, p, isCD=args.isCD, res=0)
    X_theta = build_design_matrix(theta, p, isCD=args.isCD, res=1)
    X_phi = build_design_matrix(theta, p, isCD=args.isCD, res=2)
    delta_p_corr     = X_p @ coeffs['delta_p']['beta']
    delta_theta_corr = X_theta @ coeffs['delta_theta']['beta']
    delta_phi_corr   = X_phi @ coeffs['delta_phi']['beta']

    res_p     = delta_p_orig - delta_p_corr
    res_theta = delta_theta_orig - delta_theta_corr
    res_phi   = delta_phi_orig - delta_phi_corr

    # Update figure to have 6 rows, 3 columns
    fig, axs = plt.subplots(6, 3, figsize=(18, 25))

    # # --- Row 0: original residuals heat map vs (p, theta) ---
    # for i, (res, title) in enumerate(zip([delta_p_orig, delta_theta_orig, delta_phi_orig],
    #                                     ['Δp','Δθ','Δφ'])):
    #     h = axs[0,i].hist2d(p, theta, bins=[60,60],
    #                         range=[p_range, theta_range], cmap='magma')
    #     fig.colorbar(h[3], ax=axs[0,i])
    #     axs[0,i].axhline(0, color="cyan", linestyle="--", linewidth=1.2)
    #     axs[0,i].set_xlabel('p [GeV]')
    #     axs[0,i].set_ylabel('theta [deg]')
    #     axs[0,i].set_title(f"Original {title} residuals vs (p, theta)")

    # --- Row 1: fitted grids (unchanged) ---
    fit_grids = [delta_p_grid, delta_theta_grid, delta_phi_grid]
    fit_labels = ['Δp fit', 'Δθ fit', 'Δφ fit']
    fit_keys   = ['delta_p', 'delta_theta', 'delta_phi']
    linthresh = 1e-3

    for i, (grid, label, key) in enumerate(zip(fit_grids, fit_labels, fit_keys)):
        c = axs[1,i].pcolormesh(
            P, THETA, grid, shading='auto', cmap='viridis',
            norm=SymLogNorm(linthresh=linthresh, vmin=grid.min(), vmax=grid.max())
        )
        fig.colorbar(c, ax=axs[1,i])
        axs[1,i].set_xlabel('p [GeV]')
        axs[1,i].set_ylabel('theta [deg]')
        chi2_ndof = coeffs[key]["chi2_ndof"]
        axs[1,i].set_title(f"{label}\nχ²/ndof = {chi2_ndof:.9f}")

    # --- Row 2: original residuals vs p (existing) ---
    ranges = {
        'Δp': [-0.04, 0.04] if not args.isCD else [-0.075,0.075],
        'Δθ': [-1,1] if not args.isCD else [-0.75,0.75],
        'Δφ': [-2,2] if not args.isCD else [-0.4,0.4]
    }
    for i, (res, title) in enumerate(zip([delta_p_orig, delta_theta_orig, delta_phi_orig],
                                         ['Original Δp','Original Δθ','Original Δφ'])):
        h = axs[2,i].hist2d(p, res, bins=[60,60], range=[p_range, ranges[title.split()[-1]]], cmap='magma')
        fig.colorbar(h[3], ax=axs[2,i])
        axs[2,i].axhline(0, color="cyan", linestyle="--", linewidth=1.2)
        add_profile(axs[2,i], p, res, bins=30, color="w", label="Profile mean", xrange=p_range)
        axs[2,i].set_xlabel('p [GeV]')
        axs[2,i].set_ylabel(title.split()[-1])
        axs[2,i].set_title(title)

    # --- Row 3: corrected residuals vs p (existing) ---
    for i, (res, title) in enumerate(zip([res_p, res_theta, res_phi],
                                         ['Corrected Δp','Corrected Δθ','Corrected Δφ'])):
        h = axs[3,i].hist2d(p, res, bins=[60,60], range=[p_range, ranges[title.split()[-1]]], cmap='magma')
        fig.colorbar(h[3], ax=axs[3,i])
        axs[3,i].axhline(0, color="cyan", linestyle="--", linewidth=1.2)
        add_profile(axs[3,i], p, res, bins=30, color="w", label="Profile mean", xrange=p_range)
        axs[3,i].set_xlabel('p [GeV]')
        axs[3,i].set_ylabel(title.split()[-1])
        axs[3,i].set_title(title)

    # --- Row 4: original residuals vs theta ---
    for i, (res, title) in enumerate(zip([delta_p_orig, delta_theta_orig, delta_phi_orig],
                                        ['Original Δp','Original Δθ','Original Δφ'])):
        h = axs[4,i].hist2d(theta, res, bins=[60,60],
                            range=[theta_range, ranges[title.split()[-1]]], cmap='magma')
        fig.colorbar(h[3], ax=axs[4,i])
        axs[4,i].axhline(0, color="cyan", linestyle="--", linewidth=1.2)
        add_profile(axs[4,i], theta, res, bins=30, color="w", label="Profile mean", xrange=theta_range)
        axs[4,i].set_xlabel('theta [deg]')
        axs[4,i].set_ylabel(title.split()[-1])
        axs[4,i].set_title(title + " vs theta")

    # --- Row 5: corrected residuals vs theta ---
    for i, (res, title) in enumerate(zip([res_p, res_theta, res_phi],
                                        ['Corrected Δp','Corrected Δθ','Corrected Δφ'])):
        h = axs[5,i].hist2d(theta, res, bins=[60,60],
                            range=[theta_range, ranges[title.split()[-1]]], cmap='magma')
        fig.colorbar(h[3], ax=axs[5,i])
        axs[5,i].axhline(0, color="cyan", linestyle="--", linewidth=1.2)
        add_profile(axs[5,i], theta, res, bins=30, color="w", label="Profile mean", xrange=theta_range)
        axs[5,i].set_xlabel('theta [deg]')
        axs[5,i].set_ylabel(title.split()[-1])
        axs[5,i].set_title(title + " vs theta")

    plt.tight_layout()
    if args.isCD:
        plt.savefig("pCD_corrections_with_residuals.png", dpi=300)
    else:
        plt.savefig("pFD_corrections_with_residuals.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
