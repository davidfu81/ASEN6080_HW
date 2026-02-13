function Htilde_sc = Htilde_sc_rho_rhod(t, X, Rsi)
    %%%%
    % X = [x; y; z; vx; vy; vz; ...]
    % Rsi = [Rs1, Rs2, ...] in ECEF
    %%%%
    num_station = size(Rsi,2);
    Htilde_sc = zeros([2*num_station, length(X)]);
    for k = 1:num_station
    
        % Transform to ECI
        wEtilde = tilde(Constants.WE_VEC);
        
        theta = Constants.WE*t;
        C = [cos(theta), -sin(theta), 0;
            sin(theta), cos(theta), 0;
            0, 0, 1];
        Rs = C*Rsi(:,k);
        Vs = wEtilde*Rs;
    
        % Measurements
        r = X(1:3);
        v = X(4:6);
        rho = norm(r-Rs);
        rhodot = (r-Rs)'*(v-Vs)/rho;
        
        % State partials
        Htilde_sc(2*k-1, 1:3) = (r-Rs)'/rho;
        Htilde_sc(2*k, 1:3) = (v-Vs)'/rho - rhodot/rho^2*(r-Rs)';
        Htilde_sc(2*k, 4:6) = (r-Rs)'/rho;
   
    end

end