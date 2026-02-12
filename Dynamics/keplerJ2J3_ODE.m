function dydt = keplerJ2J3_ODE(t,y, mu, J2, J3)
    % y = [r; v]
    RE = 6378;
    
    dydt = zeros(size(y));

    % dr/dt = v
    dydt(1:3) = y(4:6);

    % dv/dt = acceleration
    r = norm(y(1:3));
    z  = y(3);
    a2b = -mu/r^3 * y(1:3);
    aJ2 = 3/2*J2*RE^2*mu/r^5* ...
        [(5*z^2/r^2-1); (5*z^2/r^2-1); (5*z^2/r^2-3)].*y(1:3);
    aJ3 = 1/2*mu*J3*RE^3 * ...
        [(35*z^3/r^9-15*z/r^7)*y(1); ...
        (35*z^3/r^9-15*z/r^7)*y(2); ...
        (35*z^4/r^9-30*z^2/r^7+3/r^5)];
    dydt(4:6) = a2b + aJ2 + aJ3;
    
end
