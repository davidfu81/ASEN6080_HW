clear

sma = 10e3;
ecc = 0.001;
inc = 40;
raan = 80;
argp = 40;
tanom = 0;


X0 = oscelt2cart(sma, ecc, inc, raan, argp, tanom, Constants.MU);
dx0 = [0.1; -0.1; 0.1; -1e-4; 1e-4; 0];

Phat0 = diag([ones([3,1]); 1e-3^2*ones([3,1])]);

R = diag([0.001^2, 1e-6^2]);

station_lats = [-35.398333, 40.427222, 35.247164];
station_longs = [148.981944, 355.749444, 243.205];
Rsi = zeros([3, length(station_lats)]);
for i = 1:3
    Rsi(:,i) = lla2ecef([station_lats(i), station_longs(i), 0])/1000';
end

period = 2*pi*sqrt(sma^3/Constants.MU);
tdata = 0:10:period*15;
dyn_params = struct('mu', Constants.MU, 'J2', Constants.J2, 'J3', 0, 'RE', Constants.RE);
dyn_model = DynamicsModel_J2J3(dyn_params, 6);
meas_model = MeasurementModel_RRR(Rsi, R);

tic
Xtrue = dyn_model.integrate_eom(tdata, X0+dx0);
Ydata = meas_model.simulate_measure(tdata, Xtrue, 10);
toc

if dyn_params.J3 == 0
    save("HW4/hw4_data_J2.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", "Phat0", "R", "period");
else
    save("HW4/hw4_data_J2J3.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", "Phat0", "R", "period");
end