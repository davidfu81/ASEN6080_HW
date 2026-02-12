function [rho, rhod, el] = G_rho_rhod_el(X, Xs, R)
    if nargin < 3
        R = zeros(2,2);
    end
    RE = 6378;

    % Range and range rate
    rel_pos = X(1:3,:) - Xs(1:3,:);
    rel_vel = X(4:6,:) - Xs(4:6,:);

    rho = sqrt(sum(rel_pos.^2));
    rhod = sum(rel_pos.*rel_vel)./rho;

    % Elevation
    up_unit = Xs(1:3, :)/RE;
    el = rad2deg(asin(sum(rel_pos.*up_unit)./rho));

    % add noise to range and range rate
    if ~(all(all(R==0)))
        noise = mvnrnd(zeros([1,2]), R, size(X, 2))';
        rho = rho + noise(1,:);
        rhod = rhod + noise(2,:);
    end
    
end