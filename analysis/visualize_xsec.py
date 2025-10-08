#!/usr/bin/env python3
#type:ignore
import os
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import argparse

def extract_gr_arrays(gr):
    """Extract X, Y, EYlow/EYhigh from TGraphAsymmErrors"""
    n = gr.GetN()
    phi = np.array([gr.GetX()[i] for i in range(n)])
    y   = np.array([gr.GetY()[i] for i in range(n)])
    yerr_low  = np.array([gr.GetEYlow()[i] for i in range(n)])
    yerr_high = np.array([gr.GetEYhigh()[i] for i in range(n)])
    yerr = 0.5 * (yerr_low + yerr_high)
    return phi, y, yerr

def extract_bin_means_from_graph(gr):
    """Get <Q2>, <Xb>, <-t> from stored TParameter objects"""
    q2_mean = gr.GetListOfFunctions().FindObject("q2_mean").GetVal()
    xb_mean = gr.GetListOfFunctions().FindObject("xb_mean").GetVal()
    t_mean  = gr.GetListOfFunctions().FindObject("t_mean").GetVal()
    return q2_mean, xb_mean, t_mean

def fit_phi_graph(gr):
    """Fit TGraphAsymmErrors with A + B cos φ + C cos 2φ and return fit + covariance"""
    """Safe fit: returns pars, cov or None,None if fit fails"""
    n = gr.GetN()
    y_vals = np.array([gr.GetY()[i] for i in range(n)])
    if np.allclose(np.std(y_vals), 0):
        print(f"Skipping fit for {gr.GetName()}: zero variance")
        return None, None

    tf1 = ROOT.TF1("fit_phi", "[0] + [1]*cos(x*3.14159265/180) + [2]*cos(2*x*3.14159265/180)", 0, 360)
    mean_val = np.mean(y_vals)
    tf1.SetParameters(max(mean_val,1e-6), 0.1*abs(mean_val)+1e-6, 0.1*abs(mean_val)+1e-6)

    # Attempt fit
    try:
        fit_result = gr.Fit(tf1, "QS")  # Q=quiet, S=return fit result
        fit_ptr = gr.GetFunction("fit_phi")
        if fit_result.Status() != 0 or fit_ptr is None:
            print(f"Warning: Fit failed for graph {gr.GetName()}")
            return None, None

        pars = [fit_ptr.GetParameter(i) for i in range(3)]
        cov_matrix = np.zeros((3,3))
        cov_ptr = fit_result.GetCovarianceMatrix()
        for i in range(3):
            for j in range(3):
                cov_matrix[i,j] = cov_ptr[i,j]

        return pars, cov_matrix

    except Exception as e:
        print(f"Exception during fit of {gr.GetName()}: {e}")
        return None, None

def propagate_fit_error(phi_vals, pars, cov):
    """Propagate covariance to get fit uncertainty band"""
    phi_vals = np.array(phi_vals)
    y_fit = pars[0] + pars[1]*np.cos(np.deg2rad(phi_vals)) + pars[2]*np.cos(2*np.deg2rad(phi_vals))
    J = np.zeros((len(phi_vals),3))
    J[:,0] = 1.0
    J[:,1] = np.cos(np.deg2rad(phi_vals))
    J[:,2] = np.cos(2*np.deg2rad(phi_vals))
    y_err = np.sqrt(np.einsum("ij,jk,ik->i", J, cov, J))
    return y_fit, y_err

def plot_single_phi_bin(gr, output_dir, save=True):
    phi, y, yerr = extract_gr_arrays(gr)

    phi_width = 0.5 * (phi[1] - phi[0])        # spacing between centers
    phi_widths = np.full_like(phi, phi_width)  # same width for all bins

    y *= 1e-6
    yerr *= 1e-6

    # Extract bin means from TParameter objects
    q2_mean, xb_mean, t_mean = extract_bin_means_from_graph(gr)

    # Attempt fit
    pars, cov, chi2_ndf = None, None, None
    try:
        pars, cov = fit_phi_graph(gr)
        if pars is not None and cov is not None:
            phi_fine = np.linspace(0, 360, 400)
            y_fit, y_fit_err = propagate_fit_error(phi_fine, pars, cov)
            y_fit *= 1e-6
            y_fit_err *= 1e-6

            # Compute chi² / ndf
            # y_model, _ = propagate_fit_error(phi, pars, cov)
            # chi2 = np.sum(((y - y_model) / yerr)**2)
            # ndf = len(y) - len(pars)
            # chi2_ndf = chi2 / ndf if ndf > 0 else np.nan
    except Exception as e:
        print(f"Warning: Fit failed for {gr.GetName()} with exception: {e}")

    # Plot using object-oriented API
    fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)

    # Error bars
    ax.errorbar(phi, y, xerr=phi_widths, yerr=yerr, fmt="o", color="black",
                label="data", elinewidth=1)
    
    if pars is not None and cov is not None:
        ax.plot(phi_fine, y_fit, color="m", lw=2, label="fit")
        ax.fill_between(phi_fine, y_fit - y_fit_err, y_fit + y_fit_err,
                         color="m", alpha=0.3, label="fit error")
        # Add chi²/ndf text
        # ax.text(0.98, 0.02, f"$\\chi^2$/ndf = {chi2_ndf:.2f}",
        #          transform=ax.transAxes,
        #          verticalalignment='bottom',
        #          horizontalalignment='right',
        #          bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # Labels
    ax.set_xlabel(r"$\phi$ [deg]")
    ax.set_ylabel(r"$\left(\frac{d^2\sigma}{dtd\phi}\right)_{\gamma ^{*} p \rightarrow p' \pi ^0} \; $ [nb/GeV²]",
                fontsize=14)

    # Title
    ax.set_title(f"<Q²>={q2_mean:.3f}, <Xb>={xb_mean:.3f}, <-t>={t_mean:.3f}")

    # Bin indices text
    import re
    match = re.match(r"gr_phi_q(\d+)_xb(\d+)_t(\d+)", gr.GetName())
    if match:
        q2_idx, xb_idx, t_idx = map(int, match.groups())
        ax.text(0.02, 0.95,
                f"Bin indices:\nQ2={q2_idx}, Xb={xb_idx}, t={t_idx}",
                transform=ax.transAxes,
                verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # Limits, grid, legend
    ax.set_xlim(0, 360)
    ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.legend()

    # Save and show
    if save:
        fname = os.path.join(output_dir, f"{gr.GetName()}.png")
        fig.savefig(fname)
    plt.show()

def plot_all_phi_bins(rootfile, output_dir="xsec_phi_plots"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    f = ROOT.TFile(rootfile)
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if isinstance(obj, ROOT.TGraphAsymmErrors):
            plot_single_phi_bin(obj, output_dir)
    f.Close()
    print(f"All phi plots saved to folder '{output_dir}'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize all phi bins with fits")
    parser.add_argument("rootfile", type=str, help="phi_xsec.root file")
    parser.add_argument("--outdir", type=str, default="phi_plots", help="Folder to save figures")
    args = parser.parse_args()

    plot_all_phi_bins(args.rootfile, args.outdir)
