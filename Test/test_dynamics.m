clear
sma = 10e3;
ecc = 0.001;
inc = 40;
raan = 80;
argp = 40;
tanom = 0;

X0 = oscelt2cart(sma, ecc, inc, raan, argp, tanom, Constants.MU);

dx0 = [1; -1; 0; 0.001; 0; -0.001];

params = struct('mu', Constants.MU, 'J2', Constants.J2, 'J3', Constants.J3);

period = 2*pi*sqrt(sma^3/Constants.MU);
teval = 0:10:period*15;

model = DynamicsModel_J2J3(params, 6); 
tic
Xhist = model.integrate_eom(teval, X0+dx0);
toc

%%
tic
[Xref, Phi] = model.integrate_eomwPhi(teval, X0);
toc

XhistPhi = zeros([length(X0),length(teval)]);
for i = 1:length(teval)
    XhistPhi(:,i) = Xref(:,i) + Phi(:,:,i)*dx0;
end

%%
close all

figure
tiledlayout(6,1)
for i = 1:6
    nexttile; grid on; hold on;
    plot(teval/period, Xhist(i,:))
    % plot(teval/period, Xref(i,:))
end

figure
tiledlayout(6,1)
stop = length(teval);
for i = 1:6
    nexttile; grid on;
    plot(teval(1:stop)/period, Xhist(i,1:stop) - XhistPhi(i,1:stop))
end
sgtitle("Nonlinear Integration vs Integrated STM")