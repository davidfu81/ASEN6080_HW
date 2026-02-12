clear

% Constants
mu = 3.986e5;
J2 = 1.08264e-3;
J3 = -2.5324e-6;

load("HW2/measurement_J3.mat")
% half = floor(length(teval)/2);
% teval = teval(1:half);
% Xtrue = Xtrue(:,1:half);
% Ydata = Ydata(:,1:half,:);

% Reference state
sma = 10e3;
period = 2*pi*sqrt(sma^3/mu);
X0 = oscelt2cart(sma, 0.001, deg2rad(40), deg2rad(80), deg2rad(40), 0, mu);
Phat0 = diag([1, 1, 1, 0.001^2, 0.001^2, 0.001^2]);

R = diag([1e-3^2, 1e-6^2]);

% Run filter
station_lats = [-35.398333, 40.427222, 35.247164];
station_longs = [148.981944, 355.749444, 243.205];

tic
[Xhist, Phist] = EKF(teval, Ydata, ...
    Phat0, X0, mu, J2, Xs_all, R, 100);
toc
Yest = simulate_measure(teval, Xhist, station_lats, station_longs, 122);

%% Plots
close all

[fig_err, ax_err] = plot_filter_error(teval, Xhist, Phist, Xtrue, "EKF Post-Fit State Error");
[fig_resid, ax_resid] = plot_meas_resid(teval, Yest, Ydata, "EKF Post-Fit Measurement Residuals");
%%
% saveas(fig_err, "Figures/EKF_err_dx.png");
saveas(fig_resid, "Figures/EKF_meas_resid_J3.png");

% for i = 1:3
%     ylim(ax_err(i), [-2e6, 2e6])
%     xlim(ax_err(i), [0, 10e3])
% end
% for i = 4:6
%     ylim(ax_err(i), [-2000, 2000])
%     xlim(ax_err(i), [0, 10e3])
% end
% saveas(fig_err, "Figures/EKF_err_trans_cov.png");

for i = 1:3
    ylim(ax_err(i), [-300, 300])
    % xlim(ax_err(i), [0, 150e3])
end
for i = 4:6
    ylim(ax_err(i), [-0.3, 0.3])
    % xlim(ax_err(i), [0, 150e3])
end
saveas(fig_err, "Figures/EKF_err_ss_J3.png");

% ylim(ax_resid(1), [-5,5])
% % yticks(ax_resid(1), -5:1:5)
% ylim(ax_resid(2), [-0.005, 0.005])
% % yticks(ax_resid(2), -0.005:0.001:0.005)
% saveas(fig_resid, "Figures/EKF_meas_resid_zoom_half.png")

%%
[rms_state_comp, rms_state_3d, rms_meas] = ...
        compute_rms(Xtrue, Xhist,Ydata, Yest);

fprintf("EKF RMS Errors:\n")
fprintf("Component-wise RMS State Error: \n[ ")
for i = 1:3
    fprintf("%.3e m ", rms_state_comp(i)*1000)
end
for i = 4:6
    fprintf("%.3e m/s ", rms_state_comp(i)*1000)
end
fprintf("]\n")
fprintf("3D RMS State Error: %.3e m, %.3e m/s\n", rms_state_3d(1)*1000, rms_state_3d(2)*1000);
fprintf("RMS Measurement Errors: %.3e m, %.3e m/s\n\n", rms_meas(1)*1000, rms_meas(2)*1000);

start = find(~isnan(Ydata(1,:,1)),1);
[rms_state_comp, rms_state_3d, rms_meas] = ...
        compute_rms(Xtrue(:,start:end), Xhist(:,start:end), ...
        Ydata(:,start:end,:), Yest(:,start:end,:));

fprintf("Without First Measurement Arc:\n")
fprintf("Component-wise RMS State Error: \n[ ")
for i = 1:3
    fprintf("%.3e m ", rms_state_comp(i)*1000)
end
for i = 4:6
    fprintf("%.3e m/s ", rms_state_comp(i)*1000)
end
fprintf("]\n")
fprintf("3D RMS State Error: %.3e m, %.3e m/s\n", rms_state_3d(1)*1000, rms_state_3d(2)*1000);
fprintf("RMS Measurement Errors: %.3e m, %.3e m/s\n\n", rms_meas(1)*1000, rms_meas(2)*1000);

