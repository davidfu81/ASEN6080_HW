function dydt = keplerJ2_wPhi_ODE(t,y, mu)
    % y is length (n^2 + n) where n is 7
    dydt = zeros(size(y));

    dydt(1:6) = keplerJ2J3_ODE(t, y(1:6), mu, y(7), 0);

    dydt(7) = 0;
    
    phi = reshape(y(8:end), [7, 7]);
    A = dfdx_wJ2J3(y(1:3), y(4:6), mu, y(7), 0);
    phidot = A([1:6, 8], [1:6,8])*phi;

    dydt(8:end) = reshape(phidot, [49,1]);
end

