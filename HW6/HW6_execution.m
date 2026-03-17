clear
load("HW6/data_J2J3.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", ...
    "Phat0", "R", "period")
rms_mask = ~(all(isnan(Ydata(1,:,:)),3));

dyn_model_J2 = DynamicsModel_J2J3(struct('mu', Constants.MU, 'J2', ...
    Constants.J2, 'J3', 0, 'RE', Constants.RE), 6);
meas_model = MeasurementModel_RRR(Rsi, R);

tic
ukf = Filter_UKF(1, 2);

[Xest_ukf, Pest_ukf, dYpre_ukf, dYpost_ukf] = ukf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2, meas_model);
toc


%%
sigma = 1e-5/1e3;
dyn_model_snc = DynamicsModel_SNC(dyn_model_J2, sigma, sigma, sigma);
ukf.set_alpha(1);
tic
[Xest_ukf_snc, Pest_ukf_snc, dYpre_ukf_snc, dYpost_ukf_snc] = ukf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_snc, meas_model);
toc

%%
tic
ukf.set_alpha(1e-4);
[Xest_ukf_alpha, Pest_ukf_alpha, dYpre_ukf_alpha, dYpost_ukf_alpha] = ukf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2, meas_model);
toc

%%
tic
ekf = Filter_EKF();
[Xest_ekf, Pest_ekf, dYpre_ekf, dYpost_ekf] = ekf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_snc, meas_model, 0);
toc

%%
close all
fig_ukf_resid = plot_meas_resid(tdata, dYpost_ukf(:,:,:)*1000, "UKF Post-Fit Residuals: \sigma = 0 m/s^2, \alpha = 1");
exportgraphics(fig_ukf_resid, "Figures/HW6/ukf_resid.png")
fig_ukf_err = plot_filter_error(tdata, Xest_ukf*1000, Pest_ukf*1e6, Xtrue*1e3, "UKF State Errors: \sigma = 0 m/s^2, \alpha = 1", rms_mask);
exportgraphics(fig_ukf_err, "Figures/HW6/ukf_err.png")

fig_ukf_snc_resid = plot_meas_resid(tdata, dYpost_ukf_snc(:,:,:)*1000, "UKF Post-Fit Residuals: \sigma = 1e-5 m/s^2, \alpha = 1");
exportgraphics(fig_ukf_snc_resid, "Figures/HW6/ukf_snc_resid.png")
fig_ukf_snc_err = plot_filter_error(tdata, Xest_ukf_snc*1000, Pest_ukf_snc*1e6, Xtrue*1e3, "UKF State Errors: \sigma = 1e-5 m/s^2, \alpha = 1", rms_mask);
exportgraphics(fig_ukf_snc_err, "Figures/HW6/ukf_snc_err.png")

fig_ukf_alpha_resid = plot_meas_resid(tdata, dYpost_ukf_alpha(:,:,:)*1000, "UKF Post-Fit Residuals: \sigma = 1e-5 m/s^2, \alpha = 1e-4");
exportgraphics(fig_ukf_alpha_resid, "Figures/HW6/ukf_alpha_resid.png")
fig_ukf_alpha_err = plot_filter_error(tdata, Xest_ukf_alpha*1000, Pest_ukf_alpha*1e6, Xtrue*1e3, "UKF State Errors: \sigma = 1e-5 m/s^2, \alpha = 1e-4", rms_mask);
exportgraphics(fig_ukf_alpha_err, "Figures/HW6/ukf_alpha_err.png")

fig_ekf_resid = plot_meas_resid(tdata, dYpost_ekf(:,:,:)*1000, "EKF Post-Fit Residuals");
exportgraphics(fig_ekf_resid, "Figures/HW6/ekf_resid.png")
fig_ekf_err = plot_filter_error(tdata, Xest_ekf*1000, Pest_ekf*1e6, Xtrue*1e3, "EKF State Errors", rms_mask);
exportgraphics(fig_ekf_err, "Figures/HW6/ekf_err.png")

%%
ukf.set_alpha(1);
tic
[Xest_ukf_delta, Pest_ukf_delta, dYpre_ukf_delta, dYpost_ukf_delta] = ukf.run_filter(tdata, Ydata, X0+299*dx0, Phat0*100, ...
        dyn_model_snc, meas_model);
toc
tic
[Xest_ekf_delta, Pest_ekf_delta, dYpre_ekf_delta, dYpost_ekf_delta] = ekf.run_filter(tdata, Ydata, X0+299*dx0, Phat0*100, ...
        dyn_model_snc, meas_model, 0);
toc

fig_ukf_delta_resid = plot_meas_resid(tdata, dYpost_ukf_delta(:,:,:,end)*1000, "UKF w/ Large Initial Error Post-Fit Residuals");
exportgraphics(fig_ukf_delta_resid, "Figures/HW6/ukf_delta_resid.png")
fig_ukf_delta_err = plot_filter_error(tdata, Xest_ukf_delta*1000, Pest_ukf_delta*1e6, Xtrue*1e3, "UKF w/ Large Initial Error State Errors");
exportgraphics(fig_ukf_delta_err, "Figures/HW6/ukf_delta_err.png")

fig_ekf_delta_resid = plot_meas_resid(tdata, dYpost_ekf_delta(:,:,:,end)*1000, "EKF w/ Large Initial Error Post-Fit Residuals");
exportgraphics(fig_ekf_delta_resid, "Figures/HW6/ekf_delta_resid.png")
fig_ekf_delta_err = plot_filter_error(tdata, Xest_ekf_delta*1000, Pest_ekf_delta*1e6, Xtrue*1e3, "EKF w/ Large Initial Error State Errors");
exportgraphics(fig_ekf_delta_err, "Figures/HW6/ekf_delta_err.png")

%%
ukf.set_alpha(1)
dyn_model_J2J3 = DynamicsModel_J2J3(struct('mu', Constants.MU, 'J2', ...
    Constants.J2, 'J3', Constants.J3, 'RE', Constants.RE), 6);
tic
[Xest_ukf_J3, Pest_ukf_J3, dYpre_ukf_J3, dYpost_ukf_J3] = ukf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2J3, meas_model);
toc
tic
[Xest_ekf_J3, Pest_ekf_J3, dYpre_ekf_J3, dYpost_ekf_J3] = ekf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2J3, meas_model, 0);
toc

fig_ukf_J3_resid = plot_meas_resid(tdata, dYpost_ukf_J3(:,:,:,end)*1000, "UKF w/ J_3 Post-Fit Residuals");
exportgraphics(fig_ukf_J3_resid, "Figures/HW6/ukf_J3_resid.png")
fig_ukf_J3_err = plot_filter_error(tdata, Xest_ukf_J3*1000, Pest_ukf_J3*1e6, Xtrue*1e3, "UKF w/ J_3 State Errors");
exportgraphics(fig_ukf_J3_err, "Figures/HW6/ukf_J3_err.png")
fig_ekf_J3_resid = plot_meas_resid(tdata, dYpost_ekf_J3(:,:,:,end)*1000, "EKF w/ J_3 Post-Fit Residuals");
exportgraphics(fig_ekf_J3_resid, "Figures/HW6/ekf_J3_resid.png")
fig_ekf_J3_err = plot_filter_error(tdata, Xest_ekf_J3*1000, Pest_ekf_J3*1e6, Xtrue*1e3, "EKF w/ J_3 State Errors");
exportgraphics(fig_ekf_J3_err, "Figures/HW6/ekf_J3_err.png")