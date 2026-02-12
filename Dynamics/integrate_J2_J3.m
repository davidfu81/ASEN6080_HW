function Xsim = integrate_J2_J3(teval, X0, mu, J2, J3)
    %%%%%%
    % X = [x; y; z; vx; vy; vz;]
    %%%%%%

    opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);
    sol = ode45(@(t,y) ode_J2_J3(t,y, mu, J2, J3), ...
        [teval(1), teval(end)], X0, opts);
    Xsim = deval(sol, teval);

end