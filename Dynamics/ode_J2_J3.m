function dXdt = ode_J2_J3(t, X, params)
    % X = [r; v]
    RE = Constants.RE;
    mu = params.mu;
    J2 = params.J2;
    J3 = params.J3;
    
    dXdt = zeros(size(X));

    % dr/dt = v
    dXdt(1:3) = X(4:6);

    % dv/dt = acceleration
    r = norm(X(1:3));
    z  = X(3);
    a2b = -mu/r^3 * X(1:3);
    aJ2 = 3/2*J2*RE^2*mu/r^5* ...
        [(5*z^2/r^2-1); (5*z^2/r^2-1); (5*z^2/r^2-3)].*X(1:3);
    aJ3 = 1/2*mu*J3*RE^3 * ...
        [(35*z^3/r^9-15*z/r^7)*X(1); ...
        (35*z^3/r^9-15*z/r^7)*X(2); ...
        (35*z^4/r^9-30*z^2/r^7+3/r^5)];
    dXdt(4:6) = a2b + aJ2 + aJ3;
    
end
