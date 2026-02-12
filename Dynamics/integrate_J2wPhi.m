function [X, Phi] = integrate_J2wPhi(teval, X0, mu, J2)
    %%%%%%
    % X = [x; y; z; vx; vy; vz;]
    %%%%%%

    opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);
    n_state = length(X0);
    Phi0 = eye(n_state);
    sol = ode45(@(t,y) ode_J2wPhi(t,y, mu, J2), ...
        [teval(1), teval(end)], [X0; reshape(Phi0, [n_state^2,1])], opts);
    Xaug = deval(sol, teval);
    X = Xaug(1:n_state,:);

    Phi = zeros(n_state,n_state, length(teval));
    for i = 1:length(teval)
        Phi(:,:,i) = reshape(Xaug(n_state+1:end,i), [n_state,n_state]);
    end
end