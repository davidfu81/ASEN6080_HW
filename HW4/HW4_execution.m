clear
load("HW4/hw4_data_J2.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", ...
    "Phat0", "R", "period")
            
dyn_model_J2 = DynamicsModel_J2J3(struct('mu', Constants.MU, 'J2', ...
    Constants.J2, 'J3', 0, 'RE', Constants.RE), 6);
meas_model = MeasurementModel_RRR(Rsi, R);

ckf_smooth = Filter_CKF_Smoother();
batch = Filter_Batch();

[Xest_batch, Pest_batch, dYpre_batch, dYpost_batch] = batch.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2, meas_model, 1);

[Xest_smooth, Pest_smooth, dYpre_smooth, dYpost_smooth] = ckf_smooth.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_J2, meas_model);

%%
close all
% rms_mask = ~all(isnan(Ydata(1,:,:)), 3);
fig_batch_resid = plot_meas_resid(tdata, dYpost_batch(:,:,:,end)*1000, "Batch Post-Fit Residuals");
fig_batch_err = plot_filter_error(tdata, Xest_batch*1e3, Pest_batch*1e6, Xtrue*1e3, "Batch Error");
fig_smooth_resid = plot_meas_resid(tdata, dYpost_smooth*1000, "Smoothed CKF Post-Fit Residuals");
fig_smooth_err = plot_filter_error(tdata, Xest_smooth*1e3, Pest_smooth*1e6, Xtrue*1e3, "Smoothed CKF Error");

exportgraphics(fig_batch_resid, "Figures/HW4/batch_resid.png")
exportgraphics(fig_batch_err, "Figures/HW4/batch_err.png")
exportgraphics(fig_smooth_resid, "Figures/HW4/smooth_resid.png")
exportgraphics(fig_smooth_err, "Figures/HW4/smooth_err.png")

%%
close all
load("HW4/hw4_data_J2J3.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", ...
    "Phat0", "R", "period")
sigma = 1e-4/1e3;
dyn_model_snc = DynamicsModel_SNC(dyn_model_J2, sigma, sigma, sigma, "ECI");
[Xest_smooth_snc, Pest_smooth_snc, dYpre_smooth_snc, dYpost_smooth_snc] = ckf_smooth.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_snc, meas_model);
fig_smooth_snc_resid = plot_meas_resid(tdata, dYpost_smooth_snc*1000, "Smoothed CKF w/ SNC Post-Fit Residuals");
[fig_smooth_snc_err, ax_smooth_snc_err] = plot_filter_error(...
    tdata, Xest_smooth_snc*1e3, Pest_smooth_snc*1e6, Xtrue*1e3, "Smoothed CKF w/ SNC Errors",rms_mask);


ckf = Filter_CKF();
[Xest_ckf_snc, Pest_ckf_snc, dYpre_ckf_snc, dYpost_ckf_snc] = ckf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_snc, meas_model, 1);
fig_ckf_snc_resid = plot_meas_resid(tdata, dYpost_ckf_snc*1000, "Unsmoothed CKF w/ SNC Post-Fit Residuals");
[fig_ckf_snc_err, ax_ckf_snc_err] = plot_filter_error(...
    tdata, Xest_ckf_snc*1000, Pest_ckf_snc*1e6, Xtrue*1e3, "Unsmoothed CKF w/ SNC Errors", rms_mask);

for i = 1:3
    ylim(ax_smooth_snc_err(i), [-300, 300])
    ylim(ax_ckf_snc_err(i), [-300, 300])
end
for i = 4:6
    ylim(ax_smooth_snc_err(i), [-0.25, 0.25])
    ylim(ax_ckf_snc_err(i), [-0.25, 0.25])
end

exportgraphics(fig_smooth_snc_resid, "Figures/HW4/smooth_snc_resid.png")
exportgraphics(fig_smooth_snc_err, "Figures/HW4/smooth_snc_err.png")
exportgraphics(fig_ckf_snc_resid, "Figures/HW4/ckf_snc_resid.png")
exportgraphics(fig_ckf_snc_err, "Figures/HW4/ckf_snc_err.png")