clear

unit = 1e3;

sma = 10e3;
ecc = 0.001;
inc = 45;
raan = 80;
argp = 40;
tanom = 0;

mu = Constants.MU;
J2 = Constants.J2;

X0 = oscelt2cart(sma, ecc, inc, raan, argp, tanom, Constants.MU);
dx0 = [0.25; -0.25; 0.25; -5e-4; 5e-4; 5e-4];

Phat0 = diag([ones([3,1]); 1e-3^2*ones([3,1])]);

R = diag([0.001^2, 1e-6^2]);

station_lats = [-35.398333, 40.427222, 35.247164];
station_longs = [148.981944, 355.749444, 243.205];
Rsi = zeros([3, length(station_lats)]);
for i = 1:3
    Rsi(:,i) = lla2ecef([station_lats(i), station_longs(i), 0])/1000';
end

period = 2*pi*sqrt(sma^3/mu);
teval = 0:10:period*5;
tic
Xhist = integrate_J2_J3(teval, X0+dx0, mu, J2, 0);
Ydata = simulate_measure(teval, Xhist, Rsi, R);
toc


tic
[Xest, Pest, dYpre, dYpost] = EKF(teval, Ydata, X0, Phat0, Rsi, R, mu, J2, 0);
toc
%%
close all
plot_filter_error(teval, Xest, Pest, Xhist, "EKF Filter Error")
plot_meas_resid(teval, dYpre(:,:,:), "EKF Pre-Fit Residuals: Iteration")
plot_meas_resid(teval, dYpost(:,:,:), "Post-Fit Residuals: Iteration ")
