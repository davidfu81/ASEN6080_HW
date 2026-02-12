function dydt = phiODE(t,y, mu, J2, ref_sol)
    % y is length (n^2) where n is 7
    x = deval(ref_sol, t);
    
    phi = reshape(y, [7, 7]);
    A = dfdx_wJ2J3(x(1:3), x(4:6), mu, J2, 0);
    phidot = A([1:6, 8], [1:6,8])*phi;

    dydt = reshape(phidot, [49,1]);
end