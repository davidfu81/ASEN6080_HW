function dXaugdt = ode_J2wPhi(t,Xaug, mu, J2)
    %%%%%%
    % Xaug is length (n^2 + n) where n is 6
    % X = [x; y; z; vx; vy; vz;]
    %%%%%%

    dXaugdt = zeros(size(Xaug));

    dXaugdt(1:6) = ode_J2_J3(t, Xaug(1:6), mu, J2, 0);

    
    phi = reshape(Xaug(7:end), [6, 6]);
    A = dfdx_J2_J3(Xaug(1:6), mu, J2, 0);
    phidot = A*phi;

    dXaugdt(7:end) = reshape(phidot, [36,1]);
end