
# ---------------------------
# Canvas and output file
# ---------------------------
n_rows = 2 if args.isCD else 4
n_cols = n_bins // n_rows
canvas = ROOT.TCanvas("c_p_delta_p_vs_p_p_thetaBins", "Proton Delta p vs p in theta bins", 3000, 400 * n_rows)
canvas.Divide(n_cols, n_rows)

# ---------------------------
# Prepare reference lists for fit parameters and hists 
# ---------------------------
theta_centers = []
fitA_vals, fitB_vals, fitC_vals = [], [], []
fitA_errs, fitB_errs, fitC_errs = [], [], []
objects = []
# ---------------------------
# Loop over theta bins
# ---------------------------
for i in range(n_bins):
    tmin = theta_edges[i]
    tmax = theta_edges[i + 1]

    # Apply θ and detector cuts
    if args.isCD:
        df_cut = df_pro.Filter(f"rec.det == 2 && theta_deg > {tmin} && theta_deg < {tmax}")
        fit_formula = "[0] + [1]*x + [2]*x*x"
        y_range = (-0.2, 0.2)
    else:
        df_cut = df_pro.Filter(f"rec.det == 1 && theta_deg > {tmin} && theta_deg < {tmax}")
        fit_formula = "[0] + [1]/x"
        y_range = (-0.065, 0.2) if tmin > 35 else (-0.065, 0.065)

    # Histogram
    hname = f"h_p_deltap_p_theta_{int(tmin)}_{int(tmax)}"
    hist = df_cut.Histo2D((hname, f"Proton #Delta p vs p, {tmin:.1f} < #theta < {tmax:.1f};p_rec [GeV]; #Delta p [GeV]",
                           75, 0, 5, 75, y_range[0], y_range[1]), "rec.p", "delta_p")
    
    objects.append(hist)

    # Draw on canvas
    canvas.cd(i+1)
    pad = ROOT.gPad
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.13)
    pad.SetLogz()
    hist.Draw("COLZ")

    # Profile and fit
    prof = hist.GetValue().ProfileX()
    prof.SetLineColor(ROOT.kBlack)
    prof.SetMarkerColor(ROOT.kBlack)
    prof.SetMarkerStyle(10)
    objects.append(prof)

    # Determine fit range
    nbins = prof.GetNbinsX()
    xmin = xmax = None
    for b in range(1, nbins+1):
        if prof.GetBinEntries(b) > 0:
            x = prof.GetBinCenter(b)
            if xmin is None:
                xmin = x
            xmax = x
    # Only fit if valid range
    if xmin is not None and xmax is not None and xmax > xmin:
        fit_func = ROOT.TF1(f"f_fit_{i}", fit_formula, xmin, xmax)
        prof.Fit(fit_func, "R")

        # Only append if fit succeeded
        theta_centers.append(0.5*(tmin + tmax))
        fitA_vals.append(fit_func.GetParameter(0))
        fitA_errs.append(fit_func.GetParError(0))
        fitB_vals.append(fit_func.GetParameter(1))
        fitB_errs.append(fit_func.GetParError(1))
        if args.isCD:
            fitC_vals.append(fit_func.GetParameter(2))
            fitC_errs.append(fit_func.GetParError(2))

        objects.append(fit_func)

    hist.Draw("COLZ")
    hist.SetStats(False) 
    prof.Draw("SAME")
    fit_func.Draw("SAME")

# Save 2D histograms canvas
canvas.SaveAs("proton_delta_p_vs_p_thetaBins.png")