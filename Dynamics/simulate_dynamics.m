function Xsim = simulate_dynamics(teval, X0, mu, J2, J3)

    opts = odeset('RelTol',1e-9, 'AbsTol', 1e-9);
    sol = ode45(@(t,y) keplerJ2J3_ODE(t,y,mu,J2,J3), [teval(1), teval(end)], X0, opts);
    Xsim = deval(sol, teval);
end