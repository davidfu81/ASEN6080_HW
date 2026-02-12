clear
sma = 10e3;
ecc = 0.001;
inc = 40;
raan = 80;
argp = 40;
tanom = 0;

X0 = oscelt2cart(sma, ecc, inc, raan, argp, tanom, Constants.MU);

dx0 = [1; -1; 0; 0.001; 0; -0.001];

period = 2*pi*sqrt(sma^3/Constants.MU);
teval = 0:10:period*15;
tic
Xhist = integrate_J2_J3(teval, X0+dx0, Constants.MU, Constants.J2, 0);
toc

%%
tic
[Xref, Phi] = integrate_J2wPhi(teval, X0, Constants.MU, Constants.J2);
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