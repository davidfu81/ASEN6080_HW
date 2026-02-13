function Y = simulate_measure(teval, Xeval, Rsi, R)
    %%%%
    % Simulate nonlinear measurements for a whole trajectory
    %  X = [x; y; z; vx; vy; vz; ... ]
    % Rsi = [Rs1, Rs2, ...] in ECEF
    %%%% 
    if nargin < 3
        R = zeros(2);
    end
    
    num_stations = size(Rsi, 2);
    Y = nan(2, length(teval), num_stations);
    
    for i = 1:num_stations
        for k = 1:length(teval)
            [Y_rho_rhod, Y_el] = G_rho_rhod_el(teval(k), ...
                Xeval(:,k), Rsi(:,i), R);
            if Y_el >= 10
                Y(:,k,i) = Y_rho_rhod;
            end
        end
    end
end