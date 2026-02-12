function [Xhist, Phist] = CKF(tdata, Ydata, dxhat0, Phat0, Xref0, mu, J2, Xs, R)
    % load("HW2/measurement.mat");
    % Integrate refrence trajectory with STM
    opts = odeset('RelTol',1e-9, 'AbsTol', 1e-9);
    Xaug0 = [Xref0; J2; reshape(eye(7), [49,1])];
    sol_ref = ode45(@(t,y) keplerJ2_wPhi_ODE(t,y,mu), ...
            [tdata(1),tdata(end)], Xaug0, opts);
    Xref_aug = deval(sol_ref, tdata);
    
    % Initialize
    dxhat_prev = dxhat0;
    Phat_prev = Phat0;
    Phi_prev = eye(6);
    iprev = 1;
    Xhist = zeros([6,length(tdata)]);
    Phist = zeros([6,6,length(tdata)]);
    Xhist(:,1) = Xref0 + dxhat0;
    Phist(:,:,1) = Phat0;

    % Find observation arcs
    obs_ind = find(~all(isnan(Ydata(1,:,:)), 3));
    if obs_ind(1) == 1
        obs_ind = obs_ind(2:end);
    end

    % Iterate through observations
    for icurr = obs_ind
        
        % Time update from previous observation
        j = 1;
        while j + iprev < icurr
            Phi = reshape(Xref_aug(8:end, iprev+j), [7,7]);
            Phi = Phi(1:6, 1:6)*inv(Phi_prev);

            Xhist(:,iprev+j) = Xref_aug(1:6,iprev+j) + Phi*dxhat_prev;
            Phist(:,:,iprev+j) = Phi*Phat_prev*Phi';
            j = j+1;
        end
        Phi_curr = reshape(Xref_aug(8:end, icurr), [7,7]);
        Phi = Phi_curr(1:6, 1:6)*inv(Phi_prev);
        dxbar = Phi*dxhat_prev;
        Pbar = Phi*Phat_prev*Phi';

        % Measurement Update
        vis_ind = find(~isnan(Ydata(1,icurr,:)));     
        Y = Ydata(:,icurr, vis_ind);
        [rho, rhod] = G_rho_rhod_el(Xref_aug(1:6,icurr), Xs(:,icurr, vis_ind));
        dy = Y- [rho; rhod];
        Htilde = Htilde_sc_rho_rhod(Xref_aug(1:3,icurr), ...
            Xref_aug(4:6,icurr), Xs(1:3,icurr, vis_ind), Xs(4:6,icurr, vis_ind));
        K = Pbar*Htilde'/(Htilde*Pbar*Htilde' + R);

        dxhat = dxbar + K*(dy - Htilde*dxbar);
        Phat = (eye(6) - K*Htilde)*Pbar;

        % Record estimates and update for next iteration
        Xhist(:,icurr) = Xref_aug(1:6,icurr) + dxhat;
        Phist(:,:,icurr) = Phat;
        dxhat_prev = dxhat;
        Phat_prev = Phat;
        iprev = icurr;
        Phi_prev = Phi_curr(1:6, 1:6);
    end

    % Finish time update for rest of time if needed
    if iprev < length(tdata)
        j = 1;
        while j + iprev <= length(tdata)
            Phi = reshape(Xref_aug(8:end, iprev+j), [7,7]);
            Phi = Phi(1:6, 1:6)*inv(Phi_prev);
    
            Xhist(:,iprev+j) = Xref_aug(1:6,iprev+j) + Phi*dxhat_prev;
            Phist(:,:,iprev+j) = Phi*Phat_prev*Phi';
            j = j+1;
        end
    end
end