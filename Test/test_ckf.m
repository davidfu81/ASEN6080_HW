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
dx0 = [0.5; -0.5; 0.5; -5e-4; 5e-4; 5e-4];

Phat0 = diag([ones([3,1]); 1e-3^2*ones([3,1])]);

R = diag([0.001^2, 1e-6^2]);

station_lats = [-35.398333, 40.427222, 35.247164];
station_longs = [148.981944, 355.749444, 243.205];
Rsi = zeros([3, length(station_lats)]);
for i = 1:3
    Rsi(:,i) = lla2ecef([station_lats(i), station_longs(i), 0])/1000';
end


params = struct('mu', Constants.MU, 'J2', Constants.J2, 'J3', 0);

period = 2*pi*sqrt(sma^3/Constants.MU);
teval = 0:10:period*15;

dyn_model = DynamicsModel_J2J3(params, 6);
meas_model = MeasurementModel_RRR(Rsi, 6);

tic
Xhist = dyn_model.integrate_eom(teval, X0+dx0);
Ydata = meas_model.simulate_measure(teval, Xhist, R, 10);
toc

ckf = Filter_CKF(dyn_model, meas_model);
tic
[Xest, Pest, dYpre, dYpost] = ckf.run_filter(teval, Ydata, X0, Phat0, R, 1);
toc
%%
close all
plot_filter_error(teval, Xest, Pest, Xhist, "CKF Filter Error")
for i = 1:size(dYpre,4)
    plot_meas_resid(teval, dYpre(:,:,:,i), sprintf("Pre-Fit Residuals: Iteration %d", i))
    plot_meas_resid(teval, dYpost(:,:,:,i), sprintf("Post-Fit Residuals: Iteration %d", i))

end
%%
tic
[Xest, Pest, dYpre, dYpost] = CKF(teval, Ydata, X0, Phat0, dyn_model, meas_model, R, 1);
toc

plot_filter_error(teval, Xest, Pest, Xhist, "CKF Filter Error")
for i = 1:size(dYpre,4)
    plot_meas_resid(teval, dYpre(:,:,:,i), sprintf("Pre-Fit Residuals: Iteration %d", i))
    plot_meas_resid(teval, dYpost(:,:,:,i), sprintf("Post-Fit Residuals: Iteration %d", i))

end