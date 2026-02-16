clear
load("HW3/hw3_data.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", "Phat0", "R")
            
params = struct('mu', Constants.MU, 'J2', Constants.J2, 'J3', 0);

dynJ2J3_model = DynamicsModel_J2J3(params, 6);
sigmas = 10.^(-15:1:-2)/1e3; % km/s^2
dyn_model = DynamicsModel_SNC(dynJ2J3_model, 0, 0, 0);
meas_model = MeasurementModel_RRR(Rsi, 6);

ckf = Filter_CKF();
post_rms_meas_ckf = zeros([2,length(sigmas)]);
post_rms_3d_ckf = zeros([2,length(sigmas)]);
for i = 1:length(sigmas)
    dyn_model.set_Q(sigmas(i), sigmas(i), sigmas(i));
    [Xest, Pest, dYpre, dYpost] = ckf.run_filter(tdata, Ydata, X0, Phat0, dyn_model, meas_model, R, 1);
    
    [~, rms_state_3d, rms_meas] = compute_rms(Xtrue, Xest, dYpost);
    post_rms_meas_ckf(:,i) = rms_meas;
    post_rms_3d_ckf(:,i) = rms_state_3d;
end

close all
fig_rms_resid_ckf = plot_rms_pn(sigmas*1e3, post_rms_meas_ckf*1000, ...
    ["Range [m]", "Range-Rate [m/s"], "CKF Post-Fit RMS Measurement Residuals vs. SNC Process Noise Magnitude");
fig_rms_state_ckf = plot_rms_pn(sigmas*1e3, post_rms_3d_ckf*1000, ...
    ["Position [m]", "Velocity [m/s"], "CKF RMS 3D State Errors vs. SNC Process Noise Magnitude");


%%
sigma = 1e-4/1e3;
dyn_model.set_Q(sigma, sigma, sigma);
[Xest, Pest, dYpre, dYpost] = ckf.run_filter(tdata, Ydata, X0, Phat0, dyn_model, meas_model, R, 1);
fig_ckf_post_resid = plot_meas_resid(tdata, dYpost*1000, "CKF Post-Fit Residuals");
fig_ckf_err = plot_filter_error(tdata, Xest, Pest, Xtrue, "CKF Filter Error");

%%
ekf = Filter_EKF();
post_rms_meas_ekf = zeros([2,length(sigmas)]);
post_rms_3d_ekf = zeros([2,length(sigmas)]);
for i = 1:length(sigmas)
    dyn_model.set_Q(sigmas(i), sigmas(i), sigmas(i));
    [Xest, Pest, dYpre, dYpost] = ekf.run_filter(tdata, Ydata, X0, Phat0, dyn_model, meas_model, R, 0);
    
    [~, rms_state_3d, rms_meas] = compute_rms(Xtrue, Xest, dYpost);
    post_rms_meas_ekf(:,i) = rms_meas;
    post_rms_3d_ekf(:,i) = rms_state_3d;
end

fig_rms_resid_ekf = plot_rms_pn(sigmas*1e3, post_rms_meas_ekf*1000, ...
    ["Range [m]", "Range-Rate [m/s"], "EKF Post-Fit RMS Measurement Residuals vs. SNC Process Noise Magnitude");
fig_rms_state_ekf = plot_rms_pn(sigmas*1e3, post_rms_3d_ekf*1000, ...
    ["Position [m]", "Velocity [m/s"], "EKF RMS 3D State Errors vs. SNC Process Noise Magnitude");

%%
sigma = 1e-5/1e3;
dyn_model.set_Q(sigma, sigma, sigma);
[Xest, Pest, dYpre, dYpost] = ekf.run_filter(tdata, Ydata, X0, Phat0, dyn_model, meas_model, R, 0);
fig_ekf_post_resid = plot_meas_resid(tdata, dYpost*1000, "EKF Post-Fit Residuals");
fig_ekf_err = plot_filter_error(tdata, Xest, Pest, Xtrue, "EKF Filter Error");