classdef Filter_CKF_alt < Filter

    methods
        function obj = Filter_CKF_alt()
           
        end
        % Iteration not currently implemented
        function [Xhist, Phist, dYpre, dYpost] = run_filter(~, tdata, Ydata, Xref0, Phat0, dyn_model, meas_model, num_iter)

            n_state = length(Xref0);
            Xhist = zeros([n_state,length(tdata)]);
            Phist = zeros([n_state,n_state,length(tdata)]);
            dYpre = nan(size(Ydata));
            dYpost = nan(size(Ydata));
        
            Xstar_prev = Xref0;
            Phat_prev = Phat0;
            dxhat_prev = zeros([n_state,1]);
            iprev = 1;
            
            % Will be overwritten if there is a measurement at t = 0
            Xhist(:,1) = Xstar_prev;
            Phist(:,:,1) = Phat_prev;
        
            % Find observation arcs
            obs_ind = find(~all(isnan(Ydata(1,:,:)), 3));
        
            % Iterate through observations
            for icurr = obs_ind
                % No time update if measurement at t = 0
                if icurr == 1
                    Pbar = Phat_prev;
                    dxbar = dxhat_prev;
                    Xstar = Xstar_prev;
                    
                % Time update
                else
                    % Integrate Trajectory
                    [Xref, Phiref] = dyn_model.integrate_eomwPhi(tdata(iprev:icurr), Xstar_prev);
                    
                    % Time update from previous observation
                    j = 1;
                    while j + iprev < icurr
                        Xhist(:,iprev+j) = Xref(:,j+1);
                        Phi = Phiref(:,:,j+1)/Phiref(:,:,j);
                        Phist(:,:,iprev+j) = Phi*Phist(:,:,iprev+j-1)*Phi' + ...
                            dyn_model.process_noise_covariance(tdata(iprev+j)-tdata(iprev+j-1), Xref(:,j));
                        j = j+1;
                    end
                    dxbar = Phiref(:,:,end)*dxhat_prev;
                    Phi = Phiref(:,:,end)/Phiref(:,:,end-1);
                    Pbar = Phi*Phist(:,:,icurr-1)*Phi' + ...
                         dyn_model.process_noise_covariance(tdata(icurr)-tdata(icurr-1), Xref(:,end-1));
                    Xstar = Xref(:, end);
            
                end
                % Measurement Update
                % Identify visible stations
                stations = find(~isnan(Ydata(1,icurr,:)));
                num_station = length(stations);
    
                Y = reshape(Ydata(:,icurr,stations), [2*num_station,1]);
                Ystar = meas_model.measure(tdata(icurr), Xstar, stations);
    
                Htilde = meas_model.Htilde(tdata(icurr), Xstar, stations);
                
                Raug = kron(eye(num_station), meas_model.R);
    
                dy = Y - Ystar;
                dYpre(:,icurr,stations) = reshape(dy, [2,1,length(stations)]);
                K = Pbar*Htilde'/(Htilde*Pbar*Htilde' + Raug);
                dxhat = dxbar + K*(dy - Htilde*dxbar);
    
                M = (eye(length(Xref0)) - K*Htilde);
                Phat = M*Pbar*M' + K*Raug*K';
    
                % Post-fit residuals
                dYpost(:,icurr,stations) = reshape(dy - Htilde* dxhat, [2,1,num_station]);
                  
        
                % Record estimates and update for next iteration
                Xhist(:,icurr) = Xstar + dxhat;
                Phist(:,:,icurr) = Phat;
        
                Xstar_prev = Xstar;
                dxhat_prev = dxhat;
                Phat_prev = Phat;
                iprev = icurr;
            end
        
            % Finish time update for rest of time if needed
            if iprev < length(tdata)
                % Integrate Trajectory
                [Xref, Phiref] = dyn_model.integrate_eomwPhi(tdata(iprev:end), Xstar_prev);
                
                % Time update from previous observation
                j = 1;
                while j + iprev <= length(tdata)
                    Xhist(:,iprev+j) = Xref(:,j+1);
                    Phi = Phiref(:,:,j+1)/Phiref(:,:,j);
                    Phist(:,:,iprev+j) = Phi*Phist(:,:,iprev+j-1)*Phi' + ...
                        dyn_model.process_noise_covariance(tdata(iprev+j)-tdata(iprev+j-1), Xref(:,j));
                    j = j+1;
                end

            end
        end
    end

end