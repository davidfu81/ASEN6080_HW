function [Xhist, Phist, dYpre, dYpost] = CKF(tdata, Ydata, Xref0, Phat0, Rsi, R, mu, J2, num_iter)

    % Initialize
    Xhist = zeros([length(Xref0),length(tdata)]);
    Phist = zeros([length(Xref0),length(Xref0),length(tdata)]);
    dYpre = nan([size(Ydata),num_iter]);
    dYpost = nan([size(Ydata),num_iter]);
    
    X0 = Xref0;
    dx0 = zeros([size(Xref0)]);
    
    for iter = 1:num_iter
        % Initialize Iteration
        dxhat_prev = dx0;
        Phat_prev = Phat0;
    
        % Integrate refrence trajectory with STM
        [Xref, Phiref] = integrate_J2wPhi(tdata, X0, mu, J2);
    
        % Iterate through observations
        for i = 1:length(tdata)
            % Time Update
            if i == 1
                Phi = Phiref(1:6,1:6,i);
            else
                Phi = Phiref(1:6,1:6,i)/Phiref(1:6,1:6,i-1);
            end
            dxbar = Phi*dxhat_prev;
            Pbar = Phi*Phat_prev*Phi';
    
            % Measurement Update
            % Identify visible stations
            stations = find(~isnan(Ydata(1,i,:)));
            if ~isempty(stations)
                Y = reshape(Ydata(:,i,stations), [2*length(stations),1]);
                Ystar = G_rho_rhod_el(tdata(i), Xref(:,i), Rsi(:,stations), zeros(2));
    
                Htilde = Htilde_sc_rho_rhod(tdata(i), Xref(:,i), Rsi(:,stations));
                
                Raug = kron(eye(length(stations)), R);
    
                dy = Y - Ystar;
                dYpre(:,i,stations, iter) = reshape(dy, [2,1,length(stations)]);
                K = Pbar*Htilde'/(Htilde*Pbar*Htilde' + Raug);
                dxhat = dxbar + K*(dy - Htilde*dxbar);
    
                M = (eye(length(Xref0)) - K*Htilde);
                Phat = M*Pbar*M' + K*Raug*K';
    
                % Post-fit residuals
                dYpost(:,i,stations, iter) = reshape(dy - Htilde* dxhat, [2,1,length(stations)]);
          
            else
                dxhat = dxbar;
                Phat = Pbar;
            end
            
            % Record estimates and update for next time step
            Xhist(:,i) = Xref(:,i) + dxhat;
            Phist(:,:,i) = Phat;
            dxhat_prev = dxhat;
            Phat_prev = Phat;
        end

        dxhat0 = Phiref(:,:,end)\dxhat;
        X0 = X0 + dxhat0;
        dx0 = dx0 - dxhat0;
    end

end