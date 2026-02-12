function [Y, Xs_all] = simulate_measure(tref, Xref, ...
    station_lats, station_longs, theta0deg, R)
    
    if nargin < 6
        R = zeros(2);
    end

    Y = zeros(2, length(tref), 3);
    Xs_all = zeros(6, length(tref), 3);
    
    for i = 1:3
        Xs = station_traj(tref, station_lats(i), station_longs(i), theta0deg);
        [rho, rhod, el] = G_rho_rhod_el(Xref, Xs, R);
        Y(:,:,i) = [rho; rhod];
    
        el_mask = el >= 10;
        Y(:, ~el_mask, i) = nan;
        Xs_all(:,:,i) = Xs;
    end
end