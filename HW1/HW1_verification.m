%% Problem 1
% Part c
clear
sol_1c = jsondecode(fileread('prob1c_solution.json'));
state = sol_1c.inputs.state;
A_matrix = dfdx_wJ2J3(state.r, state.v, state.mu, state.J2, state.J3);

rel_err = abs(A_matrix - sol_1c.outputs.A_matrix.values)./abs(sol_1c.outputs.A_matrix.values);
rel_err(sol_1c.outputs.A_matrix.values == 0) = 0;
%% Problem 2
close all
clear

mu = 3.986e5;
J2 = 1.08264e-3;
J3 = -2.5324e-6;
sma = 10e3;
period = 2*pi*sqrt(sma^3/mu);
xref_true = table2array(readtable("HW1/HW1_truth.txt"));
ylabels = ["X [km]", "Y [km]", "Z [km]", "Vx [km/s]", "Vy [km/s]", "Vz [km/s]"];

% Part a
x0 = oscelt2cart(sma, 0.001, deg2rad(40), deg2rad(80), deg2rad(40), 0, mu);
tspan = [0 15*period];
teval = 0:10:tspan(2);
opts = odeset('RelTol',1e-9, 'AbsTol', 1e-9);
sol_ref = ode45(@(t,y) keplerJ2J3_ODE(t,y,mu,J2,0), tspan, x0, opts);
xref = deval(sol_ref, teval);

dx0 = [1; 0; 0; 0; 0.01; 0];
sol_pert_int = ode45(@(t,y) keplerJ2J3_ODE(t,y,mu,J2,0), tspan, x0+dx0, opts);
x_pert_int = deval(sol_pert_int, teval);
dx_pert_int = x_pert_int - xref;

figure
tiledlayout(6,1)
for i = 1:6
    nexttile; grid on; hold on;
    plot(teval/period, xref(i,:)-xref_true(:,i+1)')
    ylabel(ylabels(i))
end
xlabel("Time [Orbits]")
sgtitle("Difference Between Integrated Reference Trajectory and True Reference")

% Part b
sol_2b = jsondecode(fileread('HW1/prob2b_solution.json'));
yin = [sol_2b.inputs.X0.values; reshape(sol_2b.inputs.Phi0.values, [49,1])];
ydot = keplerJ2_wPhi_ODE(0, yin, mu);

xdot_rel_err = abs(ydot(1:7) - sol_2b.outputs.Xdot.values)./...
    abs(sol_2b.outputs.Xdot.values);
xdot_rel_err(sol_2b.outputs.Xdot.values == 0) = 0;

phidot_rel_err = abs(reshape(ydot(8:end), [7,7]) - sol_2b.outputs.Phidot.values)./...
    abs(sol_2b.outputs.Phidot.values);
phidot_rel_err(sol_2b.outputs.Phidot.values == 0) = 0;


% Part c
X0 = [x0; J2; reshape(eye(7), [49,1])];
opts = odeset('RelTol',1e-9, 'AbsTol', 1e-9);
sol_stm = ode45(@(t,y) keplerJ2_wPhi_ODE(t,y,mu), tspan, X0, opts);
Xstm = deval(sol_stm, teval);

dx_pert_stm = zeros(7, length(teval));
for i = 1:length(teval)
    phi = reshape(Xstm(8:end, i), [7,7]);
    dx_pert_stm(:,i) = phi*[dx0;0];
end

figure
tiledlayout(6,1)
for i = 1:6
    nexttile; grid on; hold on;
    plot(teval/period, dx_pert_int(i,:))
    plot(teval/period, dx_pert_stm(i,:))
    ylabel(ylabels(i))
end
xlabel("Time [Orbits]")
lgd = legend(["Integration", "STM"], 'NumColumns', 2);
lgd.Layout.Tile = 'north';
sgtitle("Perturbations from Reference Trajectory for Two Propagation Methods")

% Part d
figure
tiledlayout(6,1)
for i = 1:6
    nexttile; grid on; hold on;
    plot(teval/period, dx_pert_int(i,:) - dx_pert_stm(i,:))
    ylabel(ylabels(i))
end
xlabel("Time [Orbits]")
sgtitle("Difference in Deviation Vectors")

%% Problem 3
% Part b

sol_3b = jsondecode(fileread('HW1/prob3b_solution.json'));
H_tilde = Htilde_sc_rho_rhod(sol_3b.inputs.spacecraft_state.r, ...
    sol_3b.inputs.spacecraft_state.v,...
    sol_3b.inputs.station_state.Rs,...
    sol_3b.inputs.station_state.Vs);

rel_err = abs(H_tilde - sol_3b.outputs.Htilde.values)./abs(sol_3b.outputs.Htilde.values);
rel_err(sol_3b.outputs.Htilde.values == 0) = 0;

% Part d
sol_3d = jsondecode(fileread('HW1/prob3d_solution.json'));
H_tilde = Htilde_obs_rho_rhod(sol_3d.inputs.spacecraft_state.r, ...
    sol_3d.inputs.spacecraft_state.v,...
    sol_3d.inputs.station_state.Rs,...
    sol_3d.inputs.station_state.Vs);

rel_err = abs(H_tilde - sol_3d.outputs.Htilde.values)./abs(sol_3d.outputs.Htilde.values);
rel_err(sol_3d.outputs.Htilde.values == 0) = 0;

%% Problem 4
close all

% Part a, b
rho_rhod_el = zeros(3, length(teval), 3);
station_lats = [-35.398333, 40.427222, 35.247164];
station_longs = [148.981944, 355.749444, 243.205];

for i = 1:3
    Xs = station_traj(teval, station_lats(i), station_longs(i), 122);
    [rho, rhod, el] = G_rho_rhod_el(xref, Xs);
    rho_rhod_el(:,:,i) = [rho; rhod; el];

    el_mask = el >= 10;
    rho_rhod_el(:, ~el_mask, i) = nan;
    valid_times = teval(el_mask);
    fprintf("First and Last Measurement Times for Station %d: t = %d, %d\n",...
        i, valid_times(1), valid_times(end))
end

meas_labels = ["Range [km]", "Range Rate [km/s]", "Elevation [deg]"];

figure
tiledlayout(3,1)
for i = 1:3
    nexttile; hold on; grid on; 
    for j = 1:3
        scatter(teval/period, rho_rhod_el(i,:, j), 4)
    end
    ylabel(meas_labels(i))
end
xlabel("Time [Orbits]")
lgd = legend(["Station 1", "Station 2", "Station 3"], 'NumColumns', 3);
lgd.Layout.Tile = 'north';
sgtitle("Range and Range-Rate Measurements for Reference Trajectory")

% Part c
meas_labels = ["Range [RU]", "Doppler Shift [kHz]", "Elevation [deg]"];
ru_dop_el = rho_rhod2ru_dop(rho_rhod_el, 8.44e9);
ru_dop_el(2,:,:) = ru_dop_el(2,:,:)/1e3;
figure
tiledlayout(3,1)
for i = 1:3
    nexttile; hold on; grid on; 
    for j = 1:3
        scatter(teval/period, ru_dop_el(i,:, j), 4)
    end
    ylabel(meas_labels(i))
end
xlabel("Time [Orbits]")
lgd = legend(["Station 1", "Station 2", "Station 3"], 'NumColumns', 3);
lgd.Layout.Tile = 'north';
sgtitle("Range and Doppler Shift Measurements for Reference Trajectory")

% Part d
noisy_rhod = rho_rhod_el(2,:,:) + normrnd(0, 0.5e-6, 1, length(teval), 3);
colors = get(groot, 'defaultAxesColorOrder');
figure
hold on;
grid on;
for i = 1:3
    scatter(teval/period, rho_rhod_el(2,:, i), 4)
    scatter(teval/period, noisy_rhod(1, :, i), 4)
end
legend(["Station 1", "Noisy Station 1"...
    "Station 2", "Noisy Station 2", "Station 3", "Noisy Station 3"], ...
    'NumColumns', 3);
xlabel("Time [Orbits]")
ylabel("Range Rate [km/s]")
title("Original and Noisy Range-Rate Measurements")

figure
hold on;
grid on;
for i = 1:3
    scatter(teval/period, (noisy_rhod(1, :, i)-rho_rhod_el(2,:, i))*1e6, ...
        4, colors(i,:))
end
legend(["Station 1", "Station 2", "Station 3"], 'NumColumns', 3);
xlabel("Time [Orbits]")
ylabel("Range Rate Difference [mm/s]")
title("Noisy Range-Rate Measurements - Original Measurements")