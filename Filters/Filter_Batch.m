classdef Filter_Batch < Filter

    methods
        function obj = Filter_Batch()

        end

        function [Xhist, Phist, dYpre, dYpost] = run_filter(~, tdata, Ydata, Xref0, Phat0, dyn_model, meas_model, num_iter)
            
            % Initialize filter
            n_meas  = size(Ydata, 1);
            iter = 1;
            P0 = Phat0;
            X0 = Xref0;
            dx0 = zeros([length(X0),1]);
        
            dYpre = [];
            dYpost = [];
        
            while iter <= num_iter
                dYpre(:,:,:,iter) = nan(size(Ydata));
                dYpost(:,:,:,iter) = nan(size(Ydata));
                % Integrate refrence trajectory and state transition matrix
                [Xstar, Phi] = dyn_model.integrate_eomwPhi(tdata, X0);
        
                % Initialize iteration
                PbarR = chol(P0);
                PbarR_inv = inv(PbarR);
                Lambda = PbarR_inv*PbarR_inv';
                N = Lambda*dx0;
                
                % Iterate through observations
                for i = 1:length(tdata)
                    % Identify visible station
                    stations = find(~isnan(Ydata(1,i,:)));
                    n_stat = length(stations);
                    if isempty(stations)
                        continue
                    end
        
                    % Measurement Update
                    Y = reshape(Ydata(:,i,stations), [n_meas*n_stat,1]);
                    Ystar = meas_model.measure(tdata(i), Xstar(:,i), stations);
        
                    Htilde = meas_model.Htilde(tdata(i), Xstar(:,i), stations);
                    
                    Raug = kron(eye(n_stat), meas_model.R);
        
                    dy = Y - Ystar;
                    dYpre(:,i,stations,iter) = reshape(dy, [n_meas,1,n_stat]);
                    H = Htilde*Phi(:,:,i);
                    Lambda = Lambda + H'/Raug*H;
                    N = N + H'/Raug*dy;
                end
            
                % Solve normal equations
                LambdaR = chol(Lambda);
                LambdaR_inv = inv(LambdaR);
                P0_new = LambdaR_inv*LambdaR_inv';
                dxhat0 = P0_new*N;
        
                % Compute postfit
                for i = 1:length(tdata)
                    % Identify visible station
                    stations = find(~isnan(Ydata(1,i,:)));
                    n_stat = length(stations);
                    if isempty(stations)
                        continue
                    end
                    Htilde = meas_model.Htilde(tdata(i), Xstar(:,i), stations);
        
                    dy = reshape(dYpre(:,i,stations,iter),[n_meas*n_stat,1])- ...
                        Htilde*Phi(:,:,i)*dxhat0;
                    dYpost(:,i,stations,iter) = reshape(dy, [n_meas,1,n_stat]);
        
                end

                % Next iteration
                X0 = X0 + dxhat0;
                dx0 = dx0 - dxhat0;
                iter = iter + 1;
               
            end
        
            % Get full state and covariance history from final estimate
            [Xhist, Phi] = dyn_model.integrate_eomwPhi(tdata, X0);
            Phist = zeros([length(Xref0), length(Xref0), length(tdata)]);
            for i = 1:length(tdata)
                Phii = Phi(:,:,i);
                Phist(:,:,i) = Phii*P0_new*Phii';
            end
            
        end
    end

end