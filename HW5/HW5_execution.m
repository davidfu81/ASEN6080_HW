clear
load("HW5/data_J2.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", ...
    "Phat0", "R", "period")

dyn_model_J2 = DynamicsModel_J2J3(struct('mu', Constants.MU, 'J2', ...
    Constants.J2, 'J3', 0, 'RE', Constants.RE), 6);
meas_model = MeasurementModel_RRR(Rsi, R);

tic
srif = Filter_SRIF();
[Xest_srif, Pest_srif, dYpre_srif, dYpost_srif] = srif.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2, meas_model);
toc
tic
[Xest_srif_triu, Pest_srif_triu, dYpre_srif_triu, dYpost_srif_triu] = srif.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2, meas_model, false);
toc
tic
ckf = Filter_CKF();
[Xest_ckf, Pest_ckf, dYpre_ckf, dYpost_ckf] = ckf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2, meas_model, 1);
toc
%%
close all
fig_srif_resid = plot_meas_resid(tdata, dYpost_srif(:,:,:)*1000, "SRIF Post-Fit Residuals");
exportgraphics(fig_srif_resid, "Figures/HW5/srif_resid.png")
fig_srif_err = plot_filter_error(tdata, Xest_srif*1000, Pest_srif*1e6, Xtrue*1e3, "SRIF State Errors");
exportgraphics(fig_srif_err, "Figures/HW5/srif_err.png")

fig_srif_triu_resid = plot_meas_resid(tdata, dYpost_srif_triu(:,:,:)*1000, "SRIF w/o Upper Triangular Rbar Post-Fit Residuals");
exportgraphics(fig_srif_triu_resid, "Figures/HW5/srif_triu_resid.png")
fig_srif_triu_err = plot_filter_error(tdata, Xest_srif_triu*1000, Pest_srif_triu*1e6, Xtrue*1e3, "SRIF w/o Upper Triangular Rbar State Errors");
exportgraphics(fig_srif_triu_err, "Figures/HW5/srif_triu_err.png")

fig_ckf_resid = plot_meas_resid(tdata, dYpost_ckf(:,:,:,end)*1000, "CKF Post-Fit Residuals");
exportgraphics(fig_ckf_resid, "Figures/HW5/ckf_resid.png")
fig_ckf_err = plot_filter_error(tdata, Xest_ckf*1000, Pest_ckf*1e6, Xtrue*1e3, "CKF State Errors");
exportgraphics(fig_ckf_err, "Figures/HW5/ckf_err.png")

%%

load("HW5/data_J2J3.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", ...
    "Phat0", "R", "period")
sigma = 1e-5/1e3;
dyn_model_snc = DynamicsModel_SNC(dyn_model_J2, sigma, sigma, sigma);
tic
[Xest_srif_snc, Pest_srif_snc, dYpre_srif_snc, dYpost_srif_snc] = srif.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_snc, meas_model);
toc
tic
[Xest_ckf_snc, Pest_ckf_snc, dYpre_ckf_snc, dYpost_ckf_snc] = ckf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_snc, meas_model, 1);
toc
%%
close all
fig_srif_snc_resid = plot_meas_resid(tdata, dYpost_srif_snc(:,:,:)*1000, "SRIF w/ SNC Post-Fit Residuals");
exportgraphics(fig_srif_snc_resid, "Figures/HW5/srif_snc_resid.png")
fig_srif_snc_err = plot_filter_error(tdata, Xest_srif_snc*1000, Pest_srif_snc*1e6, Xtrue*1e3, "SRIF w/ SNC State Errors");
exportgraphics(fig_srif_snc_err, "Figures/HW5/srif_snc_err.png")

fig_ckf_snc_resid = plot_meas_resid(tdata, dYpost_ckf_snc(:,:,:,end)*1000, "CKF w/ SNC Post-Fit Residuals");
exportgraphics(fig_ckf_snc_resid, "Figures/HW5/ckf_snc_resid.png")
fig_ckf_snc_err = plot_filter_error(tdata, Xest_ckf_snc*1000, Pest_ckf_snc*1e6, Xtrue*1e3, "CKF w/ SNC State Errors");
exportgraphics(fig_ckf_snc_err, "Figures/HW5/ckf_snc_err.png")



