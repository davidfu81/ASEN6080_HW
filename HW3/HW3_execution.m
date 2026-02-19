clear
load("HW3/hw3_data.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", ...
    "Phat0", "R", "period")
            
dyn_model_J2 = DynamicsModel_J2J3(struct('mu', Constants.MU, 'J2', ...
    Constants.J2, 'J3', 0, 'RE', Constants.RE), 6);
sigmas = 10.^(-15:1:-2)/1e3; % km/s^2
dyn_model_snc = DynamicsModel_SNC(dyn_model_J2, 0, 0, 0);
meas_model = MeasurementModel_RRR(Rsi, R);

ckf = Filter_CKF();
post_rms_meas_ckf_eci = zeros([2,length(sigmas)]);
post_rms_3d_ckf_eci = zeros([2,length(sigmas)]);
tic
for i = 1:length(sigmas)
    dyn_model_snc.set_Q(sigmas(i), sigmas(i), sigmas(i));
    [Xest, Pest, dYpre, dYpost] = ckf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_snc, meas_model, 1);
    
    [~, rms_state_3d, rms_meas] = compute_rms(Xtrue, Xest, dYpost);
    post_rms_meas_ckf_eci(:,i) = rms_meas;
    post_rms_3d_ckf_eci(:,i) = rms_state_3d;
end
toc

close all
fig_rms_resid_ckf_eci = plot_rms_pn(sigmas*1e3, post_rms_meas_ckf_eci*1000, ...
    ["Range [m]", "Range-Rate [m/s]"], ...
    "CKF Post-Fit RMS Measurement Residuals vs. SNC Process Noise 1-\sigma (ECI)");
exportgraphics(fig_rms_resid_ckf_eci, "Figures/HW3/rms_resid_ckf_snc_eci.png")
fig_rms_state_ckf_eci = plot_rms_pn(sigmas*1e3, post_rms_3d_ckf_eci*1000, ...
    ["Position [m]", "Velocity [m/s]"], ...
    "CKF RMS 3D State Errors vs. SNC Process Noise 1-\sigma (ECI)");
exportgraphics(fig_rms_state_ckf_eci, "Figures/HW3/rms_state_ckf_snc_eci.png")


%%
sigma = 1e-4/1e3;
dyn_model_snc.set_Q(sigma, sigma, sigma, "ECI");
[Xest_ckf_eci, Pest_ckf_eci, dYpre_ckf_eci, dYpost_ckf_eci] = ...
    ckf.run_filter(tdata, Ydata, X0, Phat0, dyn_model_snc, meas_model, 1);

fig_ckf_snc_post_resid = plot_meas_resid(tdata, dYpost_ckf_eci*1000, ...
    sprintf("CKF Post-Fit Residuals: SNC %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3));
exportgraphics(fig_ckf_snc_post_resid, "Figures/HW3/ckf_snc_post_resid_eci.png")

fig_ckf_snc_err = plot_filter_error(tdata, Xest_ckf_eci, Pest_ckf_eci, ...
    Xtrue, sprintf("CKF Filter Error: SNC %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3));
exportgraphics(fig_ckf_snc_err, "Figures/HW3/ckf_snc_err_eci.png")
[~, rms_state_eci, rms_meas_eci] = compute_rms(Xtrue, Xest_ckf_eci, dYpost_ckf_eci)


%%
sigma = 1e-4/1e3;
sigma_x = sigma; sigma_y =  sigma/10; sigma_z = sigma;
dyn_model_snc.set_Q(sigma_x, sigma_y, sigma_z, "RIC");
[Xest_ckf_ric, Pest_ckf_ric, dYpre_ckf_ric, dYpost_ckf_ric] = ...
    ckf.run_filter(tdata, Ydata, X0, Phat0, dyn_model_snc, meas_model, 1);
fig_ckf_post_resid_ric = plot_meas_resid(tdata, dYpost_ckf_ric*1000, ...
    sprintf("CKF Post-Fit Residuals: %s = [%.1e, %.1e, %.1e] m/s^2, RIC", ...
    "\sigma", sigma_x*1e3, sigma_y*1e3, sigma_z*1e3));
exportgraphics(fig_ckf_post_resid_ric, "Figures/HW3/ckf_snc_post_resid_ric.png")
fig_ckf_err_ric = plot_filter_error(tdata, Xest_ckf_ric, Pest_ckf_ric, Xtrue, ...
    sprintf("CKF Filter Error: %s = [%.1e, %.1e, %.1e] m/s^2, RIC", "\sigma", ...
    sigma_x*1e3, sigma_y*1e3, sigma_z*1e3));
exportgraphics(fig_ckf_err_ric, "Figures/HW3/ckf_snc_err_ric.png")

[~, rms_state_ric, rms_meas_ric] = compute_rms(Xtrue, Xest_ckf_ric, dYpost_ckf_ric)

%%
ekf = Filter_EKF_slow();
post_rms_meas_ekf_snc = zeros([2,length(sigmas)]);
post_rms_3d_ekf_snc = zeros([2,length(sigmas)]);
tic
for i = 1:length(sigmas)
    dyn_model_snc.set_Q(sigmas(i), sigmas(i), sigmas(i));
    [Xest, Pest, dYpre, dYpost] = ekf.run_filter(tdata, Ydata, X0, Phat0, ...
        dyn_model_snc, meas_model, 0);
    
    [~, rms_state_3d, rms_meas] = compute_rms(Xtrue, Xest, dYpost);
    post_rms_meas_ekf_snc(:,i) = rms_meas;
    post_rms_3d_ekf_snc(:,i) = rms_state_3d;
end
toc

fig_rms_resid_ekf_eci = plot_rms_pn(sigmas*1e3, post_rms_meas_ekf_snc*1000, ...
    ["Range [m]", "Range-Rate [m/s]"], ...
    "EKF Post-Fit RMS Measurement Residuals vs. SNC Process Noise 1-\sigma (ECI)");
exportgraphics(fig_rms_resid_ekf_eci, "Figures/HW3/rms_resid_ekf_snc_eci.png")
fig_rms_state_ekf_eci = plot_rms_pn(sigmas*1e3, post_rms_3d_ekf_snc*1000, ...
    ["Position [m]", "Velocity [m/s]"], ...
    "EKF RMS 3D State Errors vs. SNC Process Noise 1-\sigma (ECI)");
exportgraphics(fig_rms_state_ekf_eci, "Figures/HW3/rms_state_ekf_snc_eci.png")

%%
sigma = 1e-4/1e3;
dyn_model_snc.set_Q(sigma, sigma, sigma, "ECI");
[Xest_ekf_eci, Pest_ekf_eci, dYpre_ekf_eci, dYpost_ekf_eci] = ...
    ekf.run_filter(tdata, Ydata, X0, Phat0, dyn_model_snc, meas_model, 0);

fig_ekf_snc_post_resid = plot_meas_resid(tdata, dYpost_ekf_eci*1000, ...
    sprintf("EKF Post-Fit Residuals: SNC %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3));
exportgraphics(fig_ekf_snc_post_resid, "Figures/HW3/ekf_snc_post_resid_eci.png")

fig_ekf_snc_err = plot_filter_error(tdata, Xest_ekf_eci, Pest_ekf_eci, Xtrue, ...
    sprintf("EKF Filter Error: SNC %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3));
exportgraphics(fig_ekf_snc_err, "Figures/HW3/ekf_snc_err_eci.png")
%% DMC
post_rms_meas_ckf_dmc = zeros([2,length(sigmas)]);
post_rms_3d_ckf_dmc = zeros([2,length(sigmas)]);
dyn_model_dmc = DynamicsModel_DMC(dyn_model_J2, 0, 0, 0, period/30);
X0_dmc = [X0; 0;0;0];
ckf = Filter_CKF_slow();

tic
for i = 1:length(sigmas)
    dyn_model_dmc.set_Q_dmc(sigmas(i), sigmas(i), sigmas(i));
    Phat0_dmc = blkdiag(Phat0, dyn_model_dmc.Q_dmc);
    [Xest, Pest, dYpre, dYpost] = ckf.run_filter(tdata, Ydata, X0_dmc, ...
        Phat0_dmc, dyn_model_dmc, meas_model, 1);
    
    [~, rms_state_3d, rms_meas] = compute_rms(Xtrue, Xest(1:6,:), dYpost);
    post_rms_meas_ckf_dmc(:,i) = rms_meas;
    post_rms_3d_ckf_dmc(:,i) = rms_state_3d;
end
toc

fig_rms_resid_ckf_dmc = plot_rms_pn(sigmas*1e3, post_rms_meas_ckf_dmc*1000, ...
    ["Range [m]", "Range-Rate [m/s]"], ...
    "CKF Post-Fit RMS Measurement Residuals vs. DMC Process Noise 1-\sigma (ECI)");
exportgraphics(fig_rms_resid_ckf_dmc, "Figures/HW3/rms_resid_ckf_dmc.png")
fig_rms_state_ckf_dmc = plot_rms_pn(sigmas*1e3, post_rms_3d_ckf_dmc*1000, ...
    ["Position [m]", "Velocity [m/s]"], ...
    "CKF RMS 3D State Errors vs. DMC Process Noise 1-\sigma (ECI)");
exportgraphics(fig_rms_state_ckf_dmc, "Figures/HW3/rms_state_ckf_dmc.png")

%%
sigma = 1e-6/1e3;
dyn_model_dmc.set_Q_dmc(sigma, sigma, sigma);
Phat0_dmc = blkdiag(Phat0, dyn_model_dmc.Q_dmc);
tic
[Xest_ckf_dmc, Pest_ckf_dmc, dYpre_ckf_dmc, dYpost_ckf_dmc] = ...
    ckf.run_filter(tdata, Ydata, X0_dmc, Phat0_dmc, dyn_model_dmc, meas_model, 1);
toc
fig_ckf_dmc_post_resid = plot_meas_resid(tdata, dYpost_ckf_dmc*1000, ...
    sprintf("CKF Post-Fit Residuals: DMC %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3));
exportgraphics(fig_ckf_dmc_post_resid, "Figures/HW3/ckf_dmc_post_resid_eci.png")

fig_ckf_dmc_err = plot_filter_error(tdata, Xest_ckf_dmc, Pest_ckf_dmc, ...
    Xtrue, sprintf("CKF Filter Error: DMC %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3));
exportgraphics(fig_ckf_dmc_err, "Figures/HW3/ckf_dmc_err_eci.png")

acc_J3_true = zeros([3,length(tdata)]);
dyn_model_J2J3 = DynamicsModel_J2J3(struct('mu', Constants.MU, 'J2', ...
    Constants.J2, 'J3', Constants.J3, 'RE', Constants.RE), 6);
for i = 1:length(tdata)
    dxdt_J2J3 = dyn_model_J2J3.eom(tdata(i), Xtrue(:,i));
    dxdt_J2 = dyn_model_J2.eom(tdata(i), Xtrue(:,i));
    acc_J3_true(:,i) = dxdt_J2J3(4:6) - dxdt_J2(4:6);
end

fig_ckf_dmc_acc = plot_dmc_acc(tdata, Xest_ckf_dmc(7:9,:), acc_J3_true, ...
    "CKF DMC Acceleration and True J3 Accelerations" );
exportgraphics(fig_ckf_dmc_acc, "Figures/HW3/ckf_dmc_acc.png")
fig_ckf_dmc_acc_err = plot_filter_error(tdata, Xest_ckf_dmc(7:9,:), ...
    Pest_ckf_dmc(7:9,7:9,:), acc_J3_true, ...
    sprintf("CKF DMC Acceleration Error: %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3),...
    ["Ax Error [m/s^2]", "Ay Error [m/s^2]", "Ax Error [m/s^2]"]);
exportgraphics(fig_ckf_dmc_acc_err, "Figures/HW3/ckf_dmc_acc_err.png")

%%
post_rms_meas_ekf_dmc = zeros([2,length(sigmas)]);
post_rms_3d_ekf_dmc = zeros([2,length(sigmas)]);

tic
for i = 1:length(sigmas)
    dyn_model_dmc.set_Q_dmc(sigmas(i), sigmas(i), sigmas(i));
    Phat0_dmc = blkdiag(Phat0, dyn_model_dmc.Q_dmc);
    [Xest, Pest, dYpre, dYpost] = ekf.run_filter(tdata, Ydata, X0_dmc, ...
        Phat0_dmc, dyn_model_dmc, meas_model, 100);
    
    [~, rms_state_3d, rms_meas] = compute_rms(Xtrue, Xest(1:6,:), dYpost);
    post_rms_meas_ekf_dmc(:,i) = rms_meas;
    post_rms_3d_ekf_dmc(:,i) = rms_state_3d;
end
toc

fig_rms_resid_ekf_dmc = plot_rms_pn(sigmas*1e3, post_rms_meas_ekf_dmc*1000, ...
    ["Range [m]", "Range-Rate [m/s]"], ...
    "EKF Post-Fit RMS Measurement Residuals vs. DMC Process Noise 1-\sigma (ECI)");
exportgraphics(fig_rms_resid_ekf_dmc, "Figures/HW3/rms_resid_ekf_dmc.png")
fig_rms_state_ekf_dmc = plot_rms_pn(sigmas*1e3, post_rms_3d_ekf_dmc*1000, ...
    ["Position [m]", "Velocity [m/s]"], ...
    "EKF RMS 3D State Errors vs. DMC Process Noise 1-\sigma (ECI)");
exportgraphics(fig_rms_state_ekf_dmc, "Figures/HW3/rms_state_ekf_dmc.png")
%%
sigma = 1e-6/1e3;
dyn_model_dmc.set_Q_dmc(sigma, sigma, sigma);
Phat0_dmc = blkdiag(Phat0, dyn_model_dmc.Q_dmc);
tic
[Xest_ekf_dmc, Pest_ekf_dmc, dYpre_ekf_dmc, dYpost_ekf_dmc] = ...
    ekf.run_filter(tdata, Ydata, X0_dmc, Phat0_dmc, dyn_model_dmc, meas_model, 100);
toc
fig_ekf_dmc_post_resid = plot_meas_resid(tdata, dYpost_ekf_dmc*1000, ...
    sprintf("EKF Post-Fit Residuals: DMC %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3));
exportgraphics(fig_ekf_dmc_post_resid, "Figures/HW3/ekf_dmc_post_resid_eci.png")

fig_ekf_dmc_err = plot_filter_error(tdata, Xest_ekf_dmc, Pest_ekf_dmc, Xtrue, ...
    sprintf("EKF Filter Error: DMC %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3));
exportgraphics(fig_ekf_dmc_err, "Figures/HW3/ekf_dmc_err_eci.png")

fig_ekf_dmc_acc = plot_dmc_acc(tdata, Xest_ekf_dmc(7:9,:), acc_J3_true, ...
    "EKF DMC Acceleration and True J3 Accelerations" );
exportgraphics(fig_ekf_dmc_acc, "Figures/HW3/ekf_dmc_acc.png")
fig_ekf_dmc_acc_err = plot_filter_error(tdata, Xest_ekf_dmc(7:9,:), ...
    Pest_ekf_dmc(7:9,7:9,:), acc_J3_true, ...
    sprintf("EKF DMC Acceleration Error: %s = %.1e m/s^2, ECI", "\sigma", sigma*1e3),...
    ["Ax Error [m/s^2]", "Ay Error [m/s^2]", "Ax Error [m/s^2]"]);
exportgraphics(fig_ekf_dmc_acc_err, "Figures/HW3/ekf_dmc_acc_err.png")

%% New tau

[~, rms_state_ckf_dmc, rms_meas_ckf_dmc] = compute_rms(Xtrue, ...
    Xest_ckf_dmc(1:6,:), dYpost_ckf_dmc)
[~, rms_state_ekf_dmc, rms_meas_ekf_dmc] = compute_rms(Xtrue, ...
    Xest_ekf_dmc(1:6,:), dYpost_ekf_dmc)


sigma = 1e-7/1e3;
dyn_model_dmc.set_Q_dmc(sigma, sigma, sigma);
dyn_model_dmc.set_tau(period/3);
[Xest_ckf_dmc_tau, Pest_ckf_dmc_tau, dYpre_ckf_dmc_tau, dYpost_ckf_dmc_tau] = ...
    ckf.run_filter(tdata, Ydata, X0_dmc, Phat0_dmc, dyn_model_dmc, meas_model, 1);
fig_ckf_dmc_post_resid_tau = plot_meas_resid(tdata, dYpost_ckf_dmc_tau*1000, ...
    sprintf("CKF Post-Fit Residuals: DMC %s = %.1e m/s^2, %s = P/3", "\sigma", ...
    sigma*1e3, "\tau"));
exportgraphics(fig_ckf_dmc_post_resid_tau, "Figures/HW3/ckf_dmc_post_resid_tau.png")

fig_ckf_dmc_err_tau = plot_filter_error(tdata, Xest_ckf_dmc_tau(1:6,:), Pest_ckf_dmc_tau, ...
    Xtrue, sprintf("CKF Filter Error: DMC %s = %.1e m/s^2, %s = P/3", "\sigma",...
    sigma*1e3, "\tau"));
exportgraphics(fig_ckf_dmc_err_tau, "Figures/HW3/ckf_dmc_err_tau.png")
[~, rms_state_ckf_tau, rms_meas_ckf_tau] = compute_rms(Xtrue, ...
    Xest_ckf_dmc_tau(1:6,:), dYpost_ckf_dmc_tau)

% EKF
[Xest_ekf_dmc_tau, Pest_ekf_dmc_tau, dYpre_ekf_dmc_tau, dYpost_ekf_dmc_tau] = ...
    ekf.run_filter(tdata, Ydata, X0_dmc, Phat0_dmc, dyn_model_dmc, meas_model, 100);
fig_ekf_dmc_post_resid_tau = plot_meas_resid(tdata, dYpost_ekf_dmc_tau*1000, ...
    sprintf("EKF Post-Fit Residuals: DMC %s = %.1e m/s^2, %s = P/3", "\sigma", ...
    sigma*1e3, "\tau"));
exportgraphics(fig_ekf_dmc_post_resid_tau, "Figures/HW3/ekf_dmc_post_resid_tau.png")

fig_ekf_dmc_err_tau = plot_filter_error(tdata, Xest_ekf_dmc_tau(1:6,:), Pest_ekf_dmc_tau, ...
    Xtrue, sprintf("EKF Filter Error: DMC %s = %.1e m/s^2, %s = P/3", "\sigma",...
    sigma*1e3, "\tau"));
exportgraphics(fig_ekf_dmc_err_tau, "Figures/HW3/ekf_dmc_err_tau.png")
[~, rms_state_ekf_tau, rms_meas_ekf_tau] = compute_rms(Xtrue, ...
    Xest_ekf_dmc_tau(1:6,:), dYpost_ekf_dmc_tau)