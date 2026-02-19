classdef Filter_EKF_slow < Filter

    methods
        function obj = Filter_EKF_slow()

        end

        function [Xhist, Phist, dYpre, dYpost] = run_filter(~, tdata, Ydata, Xref0, Phat0, dyn_model, meas_model, ckf_meas_num)
            % Initialize
            n_state = length(Xref0);
            Xhist = zeros([n_state,length(tdata)]);
            Phist = zeros([n_state,n_state,length(tdata)]);
            dYpre = nan([size(Ydata)]);
            dYpost = nan([size(Ydata)]);
                        
            % Initialize with warm start
            if ckf_meas_num > 0
                last_obs = find(~all(isnan(Ydata(1,:,:)), 3), ckf_meas_num);
                last_obs = last_obs(end);

                ckf = Filter_CKF_slow();
        
                [Xhist(:,1:last_obs), Phist(:,:,1:last_obs), dYpre(:,1:last_obs,:), dYpost(:,1:last_obs,:)] = ckf.run_filter( ...
                    tdata(1:last_obs), Ydata(:,1:last_obs,:), Xref0, Phat0, dyn_model, meas_model, 1);
        
                istart = last_obs+1;
                Xstar_prev = Xhist(:,last_obs);
                Phat_prev = Phist(:,:,last_obs);
            else
                istart = 1;
            end
        
            
            % Iterate through observations
            for i = istart:length(tdata)
                % Time Update
                if i == 1
                    Pbar = Phat0;
                    Xstar = Xref0;
                else
                    % Integrate refrence trajectory with STM
                    [Xref, Phiref] = dyn_model.integrate_eomwPhi([tdata(i-1:i)], Xstar_prev);
                    
                    Phi = Phiref(:,:,end);
                    Xstar = Xref(:,end);
                    dt = tdata(i) - tdata(i-1);
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
                    dYpre(:,i,stations) = reshape(dy, [2,1,length(stations)]);
                    K = Pbar*Htilde'/(Htilde*Pbar*Htilde' + Raug);
                    dxhat = K*dy;
        
                    M = (eye(n_state) - K*Htilde);
                    Phat = M*Pbar*M' + K*Raug*K';
        
                    % Post-fit residuals
                    dYpost(:,i,stations) = reshape(dy - Htilde* dxhat, [2,1,length(stations)]);
              
                else
                    dxhat = zeros([n_state,1]);
                    Phat = Pbar;
                end
                
                % Record estimates and update for next time step
                Xhist(:,i) = Xstar + dxhat;
                Phist(:,:,i) = Phat;
                Phat_prev = Phat;
                Xstar_prev = Xstar + dxhat;
            end
        
                
        end
    end

end