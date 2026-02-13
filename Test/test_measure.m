clear
sma = 10e3;
ecc = 0.001;
inc = 40;
raan = 120;
argp = 40;
tanom = 0;

X0 = oscelt2cart(sma, ecc, inc, raan, argp, tanom, Constants.MU);
dx0 = [1; 1; 1; -0.001; 0.001; 0.001]/1000;
alpha = [0, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000];

Ynonlin = zeros([2, length(alpha),3]);
Ylin = Ynonlin;

station_lats = [-35.398333, 40.427222, 35.247164];
station_longs = [148.981944, 355.749444, 243.205];
Rsi = zeros([3, length(station_lats)]);
for i = 1:3
    Rsi(:,i) = lla2ecef([station_lats(i), station_longs(i), 0])';
end

for k = 1:3
    Yref = G_rho_rhod_el(0, X0, Rsi(:,k), zeros(2));
    for i = 1:length(alpha)
        X = X0+alpha(i)*dx0;
        Ynonlin(:,i,k) = G_rho_rhod_el(0, X, Rsi(:,k), zeros(2));
        
        Htilde = Htilde_sc_rho_rhod(0, X0, Rsi(:,k));
        Ylin(:,i,k) = Yref + Htilde*alpha(i)*dx0;
    end
end


close all

figure
tiledlayout(2,1)
ylabels = ["Range [m]", "Range-Rate [m/s]"];
for i = 1:2
    nexttile; grid on; hold on;
    for k = 1:3
        plot(alpha, abs(Ynonlin(i,:,k)-Ylin(i,:,k))*1000)
        
    end
    ylabel(ylabels(i))
    xscale('log')
    yscale('log')
end
xlabel("Scale Factor")