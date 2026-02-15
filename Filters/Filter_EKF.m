classdef Filter_EKF < Filter

    methods
        function obj = Filter_EKF(dyn_model, meas_model)
            obj.dyn_model = dyn_model;
            obj.meas_model = meas_model;
        end

        function [Xhist, Phist, dYpre, dYpost] = run_filter(obj, tdata, Ydata, Xref0, Phat0, R, ckf_meas_num)
            if nargin < 7
                ckf_meas_num = 0;
            end
            n_state = length(Xref0);
            Xhist = zeros([n_state,length(tdata)]);
            Phist = zeros([n_state,n_state,length(tdata)]);
            dYpre = nan(size(Ydata));
            dYpost = nan(size(Ydata));
        
             % Initialize with warm start
            if ckf_meas_num > 0
                last_obs = find(~all(isnan(Ydata(1,:,:)), 3), ckf_meas_num);
                last_obs = last_obs(end);

                ckf = Filter_CKF(obj.dyn_model, obj.meas_model);
        
                [Xhist(:,1:last_obs), Phist(:,:,1:last_obs), dYpre(:,1:last_obs,:), dYpost(:,1:last_obs,:)] = ckf.run_filter( ...
                    tdata(1:last_obs), Ydata(:,1:last_obs,:), Xref0, Phat0, R, 1);
        
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
                [Xref, Phi] = obj.dyn_model.integrate_eomwPhi(tdata(iprev:icurr), Xstar_prev);
                
                % Time update from previous observation
                j = 1;
                while j + iprev < icurr
                    Xhist(:,iprev+j) = Xref(:,j+1);
                    Phist(:,:,iprev+j) = Phi(:,:,j+1)*Phat_prev*Phi(:,:,j+1)';
                    j = j+1;
                end
                Pbar = Phi(:,:,end)*Phat_prev*Phi(:,:,end)';
                Xstar = Xref(:, end);
        
                % Measurement Update
                    % Identify visible stations
                    stations = find(~isnan(Ydata(1,icurr,:)));
                    num_station = length(stations);
        
                    Y = reshape(Ydata(:,icurr,stations), [2*num_station,1]);
                    Ystar = obj.meas_model.measure(tdata(icurr), Xstar, stations, zeros(2));
        
                    Htilde = obj.meas_model.Htilde(tdata(icurr), Xstar, stations);
                    
                    Raug = kron(eye(num_station), R);
        
                    dy = Y - Ystar;
                    dYpre(:,icurr,stations) = reshape(dy, [2,1,num_station]);
                    K = Pbar*Htilde'/(Htilde*Pbar*Htilde' + Raug);
                    dxhat = K*dy;
        
                    M = (eye(n_state) - K*Htilde);
                    Phat = M*Pbar*M' + K*Raug*K';
        
                    % Post-fit residuals
                    dYpost(:,icurr,stations) = reshape(dy - Htilde* dxhat, [2,1,num_station]);
                  
        
                % Record estimates and update for next iteration
                Xhist(:,icurr) = Xstar + dxhat;
                Phist(:,:,icurr) = Phat;
        
                Xstar_prev = Xhist(:,icurr);
                Phat_prev = Phat;
                iprev = icurr;
            end
        
            % Finish time update for rest of time if needed
            if iprev < length(tdata)
                % Integrate Trajectory
                [Xref, Phi] = obj.dyn_model.integrate_eomwPhi(tdata(iprev:icurr, Xstar_prev));        
                
                % Time update from previous observation
                j = 1;
                while j + iprev <= length(tdata)
                    Xhist(:,iprev+j) = Xref(:,j);
                    Phist(:,:,iprev+j) = Phi(:,:,j)*Phat_prev*Phi(:,:,j)';
                    j = j+1;
                end
            end
        end
    end

end