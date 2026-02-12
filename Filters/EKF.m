function [Xhist, Phist] = EKF(tdata, Ydata, Phat0, Xref0, mu, J2, Xs, R, ckf_meas_num)
    if nargin < 9
        ckf_meas_num = 0;
    end

    Xhist = zeros([6,length(tdata)]);
    Phist = zeros([6,6,length(tdata)]);
    opts = odeset('RelTol',1e-9, 'AbsTol', 1e-9);

     % Initialize with warm start
    if ckf_meas_num > 0
        last_obs = find(~all(isnan(Ydata(1,:,:)), 3), ckf_meas_num);
        last_obs = last_obs(end);

        [Xhist(:,1:last_obs), Phist(:,:,1:last_obs)] = CKF( ...
            tdata(1:last_obs), Ydata(:,1:last_obs,:), zeros([6,1]), ...
            Phat0, Xref0, mu, J2, Xs, R);

        iprev = last_obs;
        Xstar_prev = Xhist(:,last_obs);
        Phat_prev = Phist(:,:,last_obs);
    else
        Xstar_prev = Xref0;
        Phat_prev = Phat0;
        iprev = 1;
        
        Xhist(:,1) = Xref0;
        Phist(:,:,1) = Phat0;
    end

    % Find observation arcs
    obs_ind = find(~all(isnan(Ydata(1,iprev+1:end,:)), 3))+iprev;

    % Iterate through observations
    for icurr = obs_ind

        % Integrate Trajectory
        Xaug0 = [Xstar_prev; J2; reshape(eye(7), [49,1])];
        sol_ref = ode45(@(t,y) keplerJ2_wPhi_ODE(t,y,mu), ...
            [tdata(iprev),tdata(icurr)], Xaug0, opts);
        
        % Time update from previous observation
        j = 1;
        Xaug = deval(sol_ref, tdata(iprev+1:icurr));
        while j + iprev < icurr
            Phi = reshape(Xaug(8:end, j), [7,7]);
            Phi = Phi(1:6, 1:6);

            Xhist(:,iprev+j) = Xaug(1:6,j);
            Phist(:,:,iprev+j) = Phi*Phat_prev*Phi';
            j = j+1;
        end
        Phi = reshape(Xaug(8:end, end), [7,7]);
        Phi = Phi(1:6, 1:6);
        Pbar = Phi*Phat_prev*Phi';
        Xstar = Xaug(1:6, end);

        % Measurement Update
        vis_ind = find(~isnan(Ydata(1,icurr,:)));     
        Y = Ydata(:,icurr, vis_ind);
        [rho, rhod] = G_rho_rhod_el(Xstar, Xs(:,icurr, vis_ind));
        dy = Y- [rho; rhod];
        Htilde = Htilde_sc_rho_rhod(Xstar(1:3), Xstar(4:6), ...
            Xs(1:3,icurr, vis_ind), Xs(4:6,icurr, vis_ind));
        K = Pbar*Htilde'/(Htilde*Pbar*Htilde' + R);

        dxhat = K*dy;
        Phat = (eye(6) - K*Htilde)*Pbar;

        % Record estimates and update for next iteration
        Xhist(:,icurr) = Xstar + dxhat;
        Phist(:,:,icurr) = Phat;

        Xstar_prev = Xhist(:,icurr);
        Phat_prev = Phat;
        iprev = icurr;
    end

    % Finish time update for rest of time if needed
    if iprev < length(tdata)
        Xaug0 = [Xstar_prev; J2; reshape(eye(7), [49,1])];
        sol_ref = ode45(@(t,y) keplerJ2_wPhi_ODE(t,y,mu), ...
                [tdata(iprev),tdata(end)], Xaug0, opts);
        j = 1;
        Xaug = deval(sol_ref, tdata(iprev+1:end));
        while j + iprev <= length(tdata)
            Phi = reshape(Xaug(8:end, j), [7,7]);
            Phi = Phi(1:6, 1:6);

            Xhist(:,iprev+j) = Xaug(1:6,j);
            Phist(:,:,iprev+j) = Phi*Phat_prev*Phi';
            j = j+1;
        end
    end
end