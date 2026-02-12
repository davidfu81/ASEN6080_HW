function [Xhist, Phist, Xhist_iter] = batch(tdata, Ydata, Phat0, Xref0, ...
    mu, J2, Xs, R, return_iter)
    
    % Initialize filter
    iter = 0;
    max_iter = 10;
    P0 = Phat0;
    X0 = Xref0;
    tspan = [tdata(1), tdata(end)];
    opts = odeset('RelTol',1e-9, 'AbsTol', 1e-9);

    if return_iter
        Xhist_iter = zeros([6,length(tdata),1]);
    end

    while iter < max_iter
        % Integrate refrence trajectory and state transition matrix
        Xstar0 = [X0; J2; reshape(eye(7), [49,1])];
        
        sol_pre = ode45(@(t,y) keplerJ2_wPhi_ODE(t,y,mu), ...
            tspan, Xstar0, opts);
        Xstar = deval(sol_pre, tdata);
    
        % Initialize iteration
        Lambda = inv(P0);
        N = zeros([6,1]);
        
        % Iterate through observations
        for i = 1:length(tdata)
            % Identify visible station
            vis_ind = find(~isnan(Ydata(1,i,:)));
            if isempty(vis_ind)
                continue
            end
            Y = Ydata(:,i, vis_ind);
        
            Htilde = Htilde_sc_rho_rhod(Xstar(1:3, i), Xstar(4:6, i), ...
                Xs(1:3,i, vis_ind), Xs(4:6,i, vis_ind));
            [rho, rhod] = G_rho_rhod_el(Xstar(1:6, i), Xs(1:6,i, vis_ind));
            dy = Y - [rho; rhod];
            Phi = reshape(Xstar(8:end, i), [7,7]);
            H = Htilde*Phi(1:6, 1:6);
    
            Lambda = Lambda + H'/R*H;
            N = N + H'/R*dy;
        end
    
        % Solve normal equations
        P0_new = inv(Lambda);
        dx0 = P0_new*N;

        X0_new = X0 + dx0;
        
        % Check convergence and update
        if norm(X0_new-X0)/norm(X0) < 1e-9
            break
        end
        X0 = X0_new;
        iter = iter + 1;

        if return_iter
             % Integrate to get full trajectory
            Xstar0 = [X0_new; J2; reshape(eye(7), [49,1])];
            sol_post = ode45(@(t,y) keplerJ2_wPhi_ODE(t,y,mu), ...
                tspan, Xstar0, opts);
            Xstar_post = deval(sol_post, tdata);
        
            Xhist = Xstar_post(1:6, :);
            Phist = zeros(6, 6, length(tdata));
            for i = 1:length(tdata)
                Phi = reshape(Xstar_post(8:end, i), [7,7]);
                Phist(:,:,i) = Phi(1:6, 1:6)*P0_new*Phi(1:6, 1:6)';
            end
            Xhist_iter(:,:,iter) = Xhist;
        end
    end

    if iter == max_iter
        fprintf("Warning: Batch filter failed to converge.")
    end
    
    if ~ return_iter
        % Integrate to get full trajectory
        Xstar0 = [X0_new; J2; reshape(eye(7), [49,1])];
        sol_post = ode45(@(t,y) keplerJ2_wPhi_ODE(t,y,mu), tspan, ...
            Xstar0, opts);
        Xstar_post = deval(sol_post, tdata);
    
        Xhist = Xstar_post(1:6, :);
        Phist = zeros(6, 6, length(tdata));
        for i = 1:length(tdata)
            Phi = reshape(Xstar_post(8:end, i), [7,7]);
            Phist(:,:,i) = Phi(1:6, 1:6)*P0_new*Phi(1:6, 1:6)';
        
        end
    end
    
end