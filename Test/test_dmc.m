close all
sma = 10e3;
ecc = 0.001;
inc = 45;
raan = 80;
argp = 40;
tanom = 0;

X0 = oscelt2cart(sma, ecc, inc, raan, argp, tanom, Constants.MU);
period = 2*pi*sqrt(sma^3/Constants.MU);
teval = 0:10:period*0.5;

params = struct('mu', Constants.MU, 'J2', Constants.J2, 'J3', Constants.J3);
sigma = 1e-5/1e3;

dyn_model_J2J3 = DynamicsModel_J2J3(params, 6);

sigma = 2e-3;
teval = 0:0.1:20;
Phat0 = blkdiag(eye(3), 1e-3*eye(3), sigma^2*eye(3));
tau = 3;
dyn_model_dmc = DynamicsModel_DMC(dyn_model_J2J3, sigma, sigma, sigma, tau);

X0 = [X0; 0; 0; 0];

dx0 = [1; -1; 0; 1e-3; 0; -1e-3; sigma*2; -2*sigma; sigma];

[Xref, Phiref] = dyn_model_dmc.integrate_eomwPhi(teval, X0);
[Xref_alt, Phiref_alt] = dyn_model_dmc.integrate_eomwPhi_alt(teval, X0);
Xnonlin = dyn_model_dmc.integrate_eom(teval, X0+dx0);
Xlin = zeros(size(Xnonlin));
for i = 1:length(teval)
    Xlin(:,i) = Phiref(:,:,i)*dx0 + Xref(:,i);
end
figure
tiledlayout(6,1)
for i = 1:6
    nexttile; grid on; hold on;
    plot(teval, Xnonlin(i,:)-Xlin(i,:), 'LineStyle', '-', 'Color', 'blue')
end

ylabels_dmc = ["Ax [m/s]", "Ay [m/s]", "Az [m/s]"];
figure
tiledlayout(3,1)
for i = 1:3
    nexttile; grid on; hold on;
    plot(teval, Xnonlin(6+i,:)-Xlin(6+i,:)*1e3, 'LineStyle', '-', 'Color', 'blue')
    ylabel(ylabels_dmc(i))
end


figure
tiledlayout(3,1)
for i = 1:3
    nexttile; grid on; hold on;
    plot(teval, Xnonlin(6+i,:)*1e3, 'LineStyle', '-', 'Color', 'blue')
    ylabel(ylabels_dmc(i))
end

Phist = zeros([size(Phat0),length(teval)]);
Phat_prev = Phat0;
for i = 1:length(teval)
    if i == 1
        Phist(:,:,i) = Phat_prev;
    else
        Phi = Phiref_alt(:,:,i)/Phiref_alt(:,:,i-1);
        if mod(i,10) == 0
            x = 1;
        end
        if mod(i,50) == 0
            x = 1;
        end
        dt = teval(i) - teval(i-1);
        Phist(:,:,i) = Phi*Phat_prev*Phi' + dyn_model_dmc.process_noise_covariance(dt, Xnonlin(:,i-1));
    end
    Phat_prev = Phist(:,:,i);
end

figure
tiledlayout(3,1)
for i = 1:3
    nexttile; grid on; hold on;
    plot(teval, Xnonlin(6+i,:)*1e3, 'LineStyle', '-', 'Color', 'blue')
    sigma = reshape(sqrt(Phist(i+6,i+6,:)), 1, []);
    plot(teval, (Xnonlin(6+i,:)-3*sigma)*1e3, 'LineStyle', '--', ...
        'Color', 'red');
    plot(teval, (Xnonlin(6+i,:)+3*sigma)*1e3, 'LineStyle', '--', ...
        'Color', 'red');
    ylabel(ylabels_dmc(i))
end