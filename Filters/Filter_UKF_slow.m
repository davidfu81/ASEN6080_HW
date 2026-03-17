classdef Filter_UKF_slow < Filter
    
    properties (Access = public)
        alpha;
        beta;
    end
    methods
        function obj = Filter_UKF_slow(alpha, beta)
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
            Xhist = zeros([length(Xref0),length(tdata)]);
            Phist = zeros([length(Xref0),length(Xref0),length(tdata)]);
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
            Xhat_prev = Xref0;
            Phat_prev = Phat0;
        
            % Iterate through observations
            for i = 1:length(tdata)

                % Time update
                if i == 1
                    Xbar = Xhat_prev;
                    Pbar = Phat_prev;
                else
                    % Compute sigma points
                    Proot = chol(Phat_prev, 'lower');
                    chi_prev = [Xhat_prev, Xhat_prev + gamma*Proot, Xhat_prev - gamma*Proot];

                    % Propagate sigma points
                    chi = zeros(size(chi_prev));
                    for j = 1:2*n_state+1
                        Xtemp = dyn_model.integrate_eom(tdata(i-1:i), chi_prev(:,j));
                        chi(:,j) = Xtemp(:,end);
                    end

                     % Compute time updated state and covariance
                    Xbar = sum(Wm.*chi,2);
                    Pbar = zeros(n_state);
                    for j = 1:n_state*2+1
                        Pbar = Pbar + Wc(j)* (chi(:,j)-Xbar)*(chi(:,j)-Xbar)';
                    end
                    Pbar = Pbar + dyn_model.process_noise_covariance(tdata(i)-tdata(i-1), Xbar);
                end
        
                % Measurement Update
                % Identify visible stations
                stations = find(~isnan(Ydata(1,i,:)));
                if ~isempty(stations)
                    % True measurement
                    num_station = length(stations);
                    Ytrue = reshape(Ydata(:,i,stations), [2*num_station,1]);

                    % Recompute sigma points
                    % Pbar = (Pbar + Pbar')/2;
                    Pbar_root = chol(Pbar, 'lower');
                    chi = [Xbar, Xbar + gamma*Pbar_root, Xbar - gamma*Pbar_root];

                    % Compute measurements from sigma poinsts
                    Ysigma = zeros([length(Ytrue), n_state*2+1]);
                    Ybar = zeros(size(Ytrue));
                    for j = 1:n_state*2+1
                        Ysigma(:,j) = meas_model.measure(tdata(i), chi(:,j), stations);
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
                    dYpre(:,i,stations) = reshape(dy, [2,1,num_station]);
                    Xhat = Xbar+K*dy;
                    Phat = Pbar - K*Pyy*K';

                    % Phat = (Phat + Phat')/2 + 1e-12*eye(n_state);
        
                    % Post-fit residuals
                    Ypost = meas_model.measure(tdata(i), Xhat, stations);
                    dYpost(:,i,stations) = reshape(Ytrue-Ypost, [2,1,num_station]);
                else
                    Xhat = Xbar;
                    Phat = Pbar;
                end    
        
                % Record estimates and update for next iteration
                Xhist(:,i) = Xhat;
                Phist(:,:,i) = Phat;
        
                Xhat_prev = Xhat;
                Phat_prev = Phat;
            end
                
        end
    end

end