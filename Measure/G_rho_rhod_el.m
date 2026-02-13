function [Y_rho_rhod, Y_el] = G_rho_rhod_el(t, X, Rsi, R)
%%%%
% Nonlinear measurement model (for one state)
% X = [x; y; z; vx; vy; vz; ...]
% Rsi = [Rs1, Rs2, ...] in ECEF
%%%%
    if nargin < 4
        R = zeros(2,2);
    end
    num_station = size(Rsi, 2);
    Y_rho_rhod = zeros([2*num_station, 1]);
    Y_el = zeros([num_station, 1]);
    for k = 1:num_station
        % Transform to ECI
        theta = Constants.WE*t;
        C = [cos(theta), -sin(theta), 0;
            sin(theta), cos(theta), 0;
            0, 0, 1];
        Rs = C*Rsi(:,k);
        Vs = cross(Constants.WE_VEC, Rs);
    
        % Range and range rate
        rel_pos = X(1:3) - Rs;
        rel_vel = X(4:6) - Vs;
    
        rho = sqrt(sum(rel_pos.^2));
        rhod = sum(rel_pos.*rel_vel)./rho;
    
        % Elevation
        up_unit = Rs/Constants.RE;
        el = rad2deg(asin(sum(rel_pos.*up_unit)./rho));
    
        % add noise to range and range rate
        if ~(all(all(R==0)))
            noise = mvnrnd(zeros([1,2]), R, size(X, 2))';
            rho = rho + noise(1,:);
            rhod = rhod + noise(2,:);
        end
        
        Y_rho_rhod(2*k-1:2*k,:) = [rho; rhod];
        Y_el(k) = el;
    end
    
end