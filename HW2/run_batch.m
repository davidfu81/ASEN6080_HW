clear

% Constants
mu = 3.986e5;
J2 = 1.08264e-3;
J3 = -2.5324e-6;

load("HW2/measurement_J3.mat")

% Reference state
sma = 10e3;
period = 2*pi*sqrt(sma^3/mu);
X0 = oscelt2cart(sma, 0.001, deg2rad(40), deg2rad(80), deg2rad(40), 0, mu);
dx0 = zeros([6,1]);
Phat0 = diag([1, 1, 1, 0.001^2, 0.001^2, 0.001^2]);

R = diag([1e-3^2, 1e-6^2]);

% Run filter
station_lats = [-35.398333, 40.427222, 35.247164];
station_longs = [148.981944, 355.749444, 243.205];

tic
half = floor(length(teval)/2);
[Xhist, Phist, Xhist_iter] = batch(teval, Ydata, ...
    Phat0, X0+dx0, mu, J2, Xs_all, R, true);
toc
Yest = simulate_measure(teval, Xhist, station_lats, station_longs, 122);

%% Plots
close all
fig_err = plot_filter_error(teval, Xhist, Phist, Xtrue, "Batch Filter Post-Fit State Error");
fig_resid = plot_meas_resid(teval, Yest, Ydata, "Batch Filter Post-Fit Measurement Residuals");
%%
if ~exist('Figures', 'dir')
    mkdir('Figures')
end
saveas(fig_err, "Figures/batch_err_J3.png");
saveas(fig_resid, "Figures/batch_meas_resid_J3.png");
%%
for iter = 1:size(Xhist_iter, 3)
    Yest = simulate_measure(teval, Xhist_iter(:,:,iter), station_lats, station_longs, 122);

    [rms_state_comp, rms_state_3d, rms_meas] = ...
        compute_rms(Xtrue, Xhist_iter(:,:,iter), Ydata, Yest);

    fprintf("Iteration %d:\n", iter)
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
end

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
fprintf("RMS Measurement Errors: %.3e m, %.3e m/s\n", rms_meas(1)*1000, rms_meas(2)*1000);