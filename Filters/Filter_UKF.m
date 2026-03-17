classdef Filter_UKF < Filter
    
    properties (Access = public)
        alpha;
        beta;
    end
    methods
        function obj = Filter_UKF(alpha, beta)
            obj.alpha = alpha;
            obj.beta = beta;
        end

        function obj = set_alpha(obj, alpha)
            obj.alpha = alpha;
        end

        function obj = set_beta(obj, beta)
            obj.beta = beta;
        end

        function [Xhist, Phist, dYpre, dYpost] = run_filter(obj, tdata, Ydata, Xref0, Phat0, dyn_model, meas_model)
            % Initialize
            Xhist = nan([length(Xref0),length(tdata)]);
            Phist = nan([length(Xref0),length(Xref0),length(tdata)]);
            dYpre = nan([size(Ydata)]);
            dYpost = nan([size(Ydata)]);

            % Unscented parameters and weights
            n_state = length(Xref0);
            kappa =  3-n_state;
            lambda = obj.alpha^2*(n_state+kappa)-n_state;
            gamma = sqrt(n_state+lambda);
            
            Wm = ones([1, 2*n_state+1])/(2*(n_state+lambda));
            Wc = Wm;
            Wm(1) = lambda/(n_state+lambda);
            Wc(1) = lambda/(n_state+lambda) + (1-obj.alpha^2 + obj.beta);
            
            % Initialize
            Xhist(:,1) = Xref0;
            Phist(:,:,1) = Phat0;
            Xhat_prev = Xref0;
            Phat_prev = Phat0;
            iprev = 1;

            % Find observation arcs
            obs_ind = find(~all(isnan(Ydata(1,1:end,:)), 3));
        
            % Iterate through observations
            for icurr = obs_ind

                % Time update
                if icurr == 1
                    Xbar = Xhat_prev;
                    Pbar = Phat_prev;
                    
                else
                    % Compute sigma points
                    try
                        Proot = chol(Phat_prev, 'lower');
                    catch
                        [V,D] = eig(Phat_prev);
                        D(D<0) = 1e-12;
                        Phat_prev = V*D/V;
                        Phist(:,:,iprev) = Phat_prev;
                        Proot = chol(Phat_prev, 'lower');
                    end
                    chi_prev = [Xhat_prev, Xhat_prev + gamma*Proot, Xhat_prev - gamma*Proot];

                    % Propagate sigma points
                    chi_hist = zeros([size(chi_prev), icurr-iprev]);
                    for j = 1:2*n_state+1
                        Xtemp = dyn_model.integrate_eom(tdata(iprev:icurr), chi_prev(:,j));
                        chi_hist(:,j,:) = Xtemp(:,2:end);
                    end

                    % Propagate time updated state and covariance
                    k = 1;
                    while k + iprev <= icurr
                        Xbar = sum(Wm.*chi_hist(:,:,k),2);
                        Pbar = zeros(n_state);
                        for j = 1:n_state*2+1
                            Pbar = Pbar + Wc(j)* (chi_hist(:,j,k)-Xbar)*(chi_hist(:,j,k)-Xbar)';
                        end
                        
                        % Add process noise if the gap between measurements is
                        % small
                        if tdata(iprev+k) - tdata(iprev) < 100
                            Pbar = Pbar + dyn_model.process_noise_covariance(tdata(iprev+k)-tdata(iprev), Xbar);
                        end

                        % Record histories
                        Xhist(:,iprev+k) = Xbar;
                        Phist(:,:,iprev+k) = Pbar;
                        k = k + 1;
                    end

                end
        
                % Measurement Update
                % Identify visible stations
                stations = find(~isnan(Ydata(1,icurr,:)));

                % True measurement
                num_station = length(stations);
                Ytrue = reshape(Ydata(:,icurr,stations), [2*num_station,1]);

                % Recompute sigma points
                Pbar = (Pbar+Pbar')/2;% + 1e-12*eye(n_state);
                try
                    Pbar_root = chol(Pbar, 'lower');
                catch
                    [V, D] = eig(Pbar);
                    D(D<0) = 1e-12;
                    Pbar = V*D/V;
                    Pbar_root = chol(Pbar, 'lower');
                end
                chi = [Xbar, Xbar + gamma*Pbar_root, Xbar - gamma*Pbar_root];

                % Compute measurements from sigma poinsts
                Ysigma = zeros([length(Ytrue), n_state*2+1]);
                Ybar = zeros(size(Ytrue));
                for j = 1:n_state*2+1
                    Ysigma(:,j) = meas_model.measure(tdata(icurr), chi(:,j), stations);
                    Ybar = Ybar + Wm(j) * Ysigma(:,j);
                end
                
                % Compute innovation and cross covariances
                Pyy = kron(eye(num_station), meas_model.R);
                Pxy = zeros([n_state, length(Ytrue)]);
                for j = 1:n_state*2+1
                    Pyy = Pyy + Wc(j)*(Ysigma(:,j)-Ybar)*(Ysigma(:,j)-Ybar)';
                    Pxy = Pxy + Wc(j)*(chi(:,j) - Xbar)*(Ysigma(:,j)-Ybar)';
                end

                % Kalman gain and measurement update;
                K = Pxy/Pyy;
                dy = Ytrue - Ybar;
                dYpre(:,icurr,stations) = reshape(dy, [2,1,num_station]);
                Xhat = Xbar+K*dy;
                Phat = Pbar - K*Pyy*K';

                Phat = (Phat + Phat')/2;
    
                % Post-fit residuals
                Ypost = meas_model.measure(tdata(icurr), Xhat, stations);
                dYpost(:,icurr,stations) = reshape(Ytrue-Ypost, [2,1,num_station]);
  
        
                % Record estimates and update for next iteration
                Xhist(:,icurr) = Xhat;
                Phist(:,:,icurr) = Phat;
        
                Xhat_prev = Xhat;
                Phat_prev = Phat;
                iprev = icurr;
            end
                
        end
    end

end