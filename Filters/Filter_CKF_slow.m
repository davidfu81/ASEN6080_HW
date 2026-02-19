classdef Filter_CKF_slow < Filter

    methods
        function obj = Filter_CKF_slow()

        end

        function [Xhist, Phist, dYpre, dYpost] = run_filter(~, tdata, Ydata, Xref0, Phat0, dyn_model, meas_model, num_iter)
            % Initialize
            n_state = length(Xref0);
            Xhist = zeros([n_state,length(tdata)]);
            Phist = zeros([n_state,n_state,length(tdata)]);
            dYpre = nan([size(Ydata),num_iter]);
            dYpost = nan([size(Ydata),num_iter]);
            
            X0 = Xref0;
            dx0 = zeros([n_state, 1]);
            
            for iter = 1:num_iter                
            
                % Iterate through observations
                for i = 1:length(tdata)
                    % Time Update
                    if i == 1
                        dxbar = dx0;
                        Pbar = Phat0;
                        Xstar = X0;
                        Phi_tot = eye(n_state);
                    else
                        % Integrate refrence trajectory with STM
                        [Xref, Phiref] = dyn_model.integrate_eomwPhi([tdata(i-1:i)], Xstar_prev);
                        
                        Phi = Phiref(:,:,end);
                        Phi_tot = Phi*Phi_tot;
                        Xstar = Xref(:,end);
                        dt = tdata(i) - tdata(i-1);
                        dxbar = Phi*dxhat_prev;
                        Pbar = Phi*Phat_prev*Phi' + dyn_model.process_noise_covariance(dt, Xstar_prev);
                    end
                    
            
                    % Measurement Update
                    % Identify visible stations
                    stations = find(~isnan(Ydata(1,i,:)));
                    if ~isempty(stations)
                        Y = reshape(Ydata(:,i,stations), [2*length(stations),1]);
                        Ystar = meas_model.measure(tdata(i), Xstar, stations);
            
                        Htilde = meas_model.Htilde(tdata(i), Xstar, stations);
                        
                        Raug = kron(eye(length(stations)), meas_model.R);
            
                        dy = Y - Ystar;
                        dYpre(:,i,stations, iter) = reshape(dy, [2,1,length(stations)]);
                        K = Pbar*Htilde'/(Htilde*Pbar*Htilde' + Raug);
                        dxhat = dxbar + K*(dy - Htilde*dxbar);
            
                        M = (eye(n_state) - K*Htilde);
                        Phat = M*Pbar*M' + K*Raug*K';
            
                        % Post-fit residuals
                        dYpost(:,i,stations, iter) = reshape(dy - Htilde* dxhat, [2,1,length(stations)]);
                  
                    else
                        dxhat = dxbar;
                        Phat = Pbar;
                    end
                    
                    % Record estimates and update for next time step
                    Xhist(:,i) = Xstar + dxhat;
                    Phist(:,:,i) = Phat;
                    dxhat_prev = dxhat;
                    Phat_prev = Phat;
                    Xstar_prev = Xstar;
                end
        
                dxhat0 = [Phi_tot(1:6,1:6)\dxhat(1:6); zeros([n_state-6,1])];
                X0 = X0 + dxhat0;
                dx0 = dx0 - dxhat0;
            end
        end
    end

end