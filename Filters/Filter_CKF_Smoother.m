classdef Filter_CKF_Smoother < Filter

    methods
        function obj = Filter_CKF_Smoother()

        end

        function [Xhist, Phist, dYpre, dYpost] = run_filter(~, tdata, Ydata, Xref0, Phat0, dyn_model, meas_model)
            % Initialize
            Xhist = zeros([length(Xref0),length(tdata)]);
            Phist = zeros([length(Xref0),length(Xref0),length(tdata)]);
            dYpre = nan(size(Ydata));
            dYpost = nan(size(Ydata));
            
            % Record histories for smoothing
            Phat_hist = Phist;
            Pbar_hist = Phist;
            dxhat_hist = Xhist;
            Htilde_hist = cell([length(tdata),1]);
            stations_hist = cell([length(tdata),1]);
            
            X0 = Xref0;
            dx0 = zeros([size(Xref0)]);
            
            % Initialize
            dxhat_prev = dx0;
            Phat_prev = Phat0;
        
            % Integrate refrence trajectory with STM
            [Xref, Phiref] = dyn_model.integrate_eomwPhi(tdata, X0);
        
            % Iterate through observations
            for i = 1:length(tdata)
                % Time Update
                if i == 1
                    dxbar = dxhat_prev;
                    Pbar = Phat_prev;
                else
                    Phi = Phiref(:,:,i)/Phiref(:,:,i-1);
                    dt = tdata(i) - tdata(i-1);
                    dxbar = Phi*dxhat_prev;
                    Pbar = Phi*Phat_prev*Phi' + dyn_model.process_noise_covariance(dt, Xref(:,i-1));
                end
                
        
                % Measurement Update
                % Identify visible stations
                stations = find(~isnan(Ydata(1,i,:)));
                if ~isempty(stations)
                    Y = reshape(Ydata(:,i,stations), [2*length(stations),1]);
                    Ystar = meas_model.measure(tdata(i), Xref(:,i), stations);
        
                    Htilde = meas_model.Htilde(tdata(i), Xref(:,i), stations);
                    
                    Raug = kron(eye(length(stations)), meas_model.R);
        
                    dy = Y - Ystar;
                    dYpre(:,i,stations) = reshape(dy, [2,1,length(stations)]);
                    K = Pbar*Htilde'/(Htilde*Pbar*Htilde' + Raug);
                    dxhat = dxbar + K*(dy - Htilde*dxbar);
        
                    M = (eye(length(Xref0)) - K*Htilde);
                    Phat = M*Pbar*M' + K*Raug*K';
        
                    % Post-fit residuals
                    % dYpost(:,i,stations) = reshape(dy - Htilde* dxhat, [2,1,length(stations)]);
              
                else
                    Htilde = [];
                    dxhat = dxbar;
                    Phat = Pbar;
                end
                
                % Record and update for next time step
                dxhat_hist(:,i) = dxhat;
                Phat_hist(:,:,i) = Phat;
                Pbar_hist(:,:,i) = Pbar;
                Htilde_hist{i} = Htilde;
                stations_hist{i} = stations;
                dxhat_prev = dxhat;
                Phat_prev = Phat;
            end
            
            dx_smooth_hist = zeros(size(dxhat_hist));
            dx_smooth_hist(:,end) = dxhat_hist(:,end);
            Phist(:,:,end) = Phat_hist(:,:,end);
            % Backwards Smoothing
            for k = length(tdata)-1:-1:1
                dx_kk = dxhat_hist(:,k);
                dx_lk1 = dx_smooth_hist(:,k+1);
                P_kk = Phat_hist(:,:,k);
                P_kk1 = Pbar_hist(:,:,k+1);
                P_lk1 = Phist(:,:,k+1);
                Phi = Phiref(:,:,k+1)/Phiref(:,:,k);
                
                dt = tdata(k) - tdata(k+1);
                Xk = Xref(:,k) + dx_kk;
                Sk = P_kk*Phi'/ (Phi*P_kk*Phi' + dyn_model.process_noise_covariance(dt, Xk));
                dx_lk = dx_kk + Sk*(dx_lk1-Phi*dx_kk);

                dx_smooth_hist(:,k) = dx_lk;
                Phist(:,:,k) = P_kk + Sk*(P_lk1 - P_kk1)*Sk';

                % Re-compute post-fit residuals
                stations = stations_hist{k};
                if ~isempty(stations)
                    Htilde = Htilde_hist{k};
                    dYpost(:,k,stations) = dYpre(:,k,stations) - reshape(Htilde* dx_lk, [2,1,length(stations)]);
                end
            end
        
            Xhist = Xref + dx_smooth_hist;
        end
    end

end