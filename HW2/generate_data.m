%% Generate Measurement Data
clear
% Constants
mu = 3.986e5;
J2 = 1.08264e-3;
J3 = -2.5324e-6;

% Trajectory
sma = 10e3;
period = 2*pi*sqrt(sma^3/mu);
ylabels = ["X [km]", "Y [km]", "Z [km]", "Vx [km/s]", "Vy [km/s]", "Vz [km/s]"];
X0 = oscelt2cart(sma, 0.001, deg2rad(40), deg2rad(80), deg2rad(40), 0, mu);
dx0 = [0.1; -0.1; 0; -0.0001; 0.0001; 0];
teval = 0:10:15*period;
Xtrue = simulate_dynamics(teval, X0+dx0, mu, J2, J3);


% Measurements
R = diag([1e-3^2, 1e-6^2]);
station_lats = [-35.398333, 40.427222, 35.247164];
station_longs = [148.981944, 355.749444, 243.205];

[Ydata, Xs_all] = simulate_measure(teval, Xtrue, station_lats, station_longs, 122, R);

save("HW2/measurement_J3.mat", "teval", "Ydata", "Xs_all", "Xtrue")

%% 
clear
close all
load("HW2/measurement_J3.mat")
Xtrue_J3 = Xtrue;
Ydata_J3 = Ydata;
load("HW2/measurement.mat")

fig_meas = plot_meas_resid(teval, Ydata_J3, Ydata, "Measurement Differences: J3 - No J3");

fig_state = figure;
tl = tiledlayout(6,1);
tl.TileSpacing = 'compact';
tl.Padding = 'loose';
ylabels = ["X Error [m]", "Y Error [m]", "Z Error [m]", ...
    "Vx Error [m/s]", "Vy Error [m/s]", "Vz Error [m/s]"];
for i = 1:6
    nexttile
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    hold on; grid on;
    errorsi = Xtrue_J3(i,:) - Xtrue(i,:);
    plot(teval, errorsi*1000, 'LineStyle', '-', ...
        'Color', 'blue')
    ylabel(ylabels(i))
end
xlabel("Time [s]")
sgtitle("State Differences: J3 - No J3");

saveas(fig_meas, "Figures/meas_compare.png")
saveas(fig_state, "Figures/state_compare.png")

