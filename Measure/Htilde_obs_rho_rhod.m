function Htilde = Htilde_obs_rho_rhod(r, v, Rs, Vs)
    % Measurements
    rho = norm(r-Rs);
    rhodot = (r-Rs)'*(v-Vs)/rho;

    % Partials
    Htilde = zeros(2,3);
    Htilde(1, 1:3) = -(r-Rs)'/rho;
    Htilde(2,1:3) = -(v-Vs)'/rho + rhodot/rho^2*(r-Rs)';
end