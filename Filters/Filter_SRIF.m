classdef Filter_SRIF < Filter

    methods
        function obj = Filter_SRIF()

        end

        function [Xhist, Phist, dYpre, dYpost] = run_filter(~, tdata, Ydata, Xref0, Phat0, dyn_model, meas_model, triu_rbar)
            if nargin < 8
                triu_rbar = true;
            end
            % Initialize
            n_state = length(Xref0);
            Xhist = zeros([n_state,length(tdata)]);
            Phist = zeros([n_state,n_state,length(tdata)]);
            dYpre = nan([size(Ydata)]);
            dYpost = nan([size(Ydata)]);
            
            X0 = Xref0;
            dx0 = zeros([n_state, 1]);
            Phat0_root = chol(Phat0);
            Rbar0 = inv(Phat0_root);

            % Integrate refrence trajectory with STM
            [Xref, Phiref] = dyn_model.integrate_eomwPhi(tdata, X0);
                        
            % Iterate through observations
            for i = 1:length(tdata)
                % Time Update
                if i == 1
                    dxbar = dx0;
                    Rbar = householder(Rbar0);
                    bbar = Rbar*dxbar;
                    Xstar = X0;
                else
                    Phi = Phiref(:,:,i)/Phiref(:,:,i-1);
                    Xstar = Xref(:,i);
                    dxbar = Phi*dxhat_prev;
                    
                    if isa(dyn_model, "DynamicsModel_SNC")
                        Q_root = chol(dyn_model.process_noise(Xref(:,i-1)+dxhat_prev));
                        Ru = inv(Q_root);
                        Rtilde = Rhat_prev/Phi;
                        
                        dt = tdata(i)-tdata(i-1);
                        n_noise = size(Ru,1);
                        A = [Ru, zeros(n_noise, n_state), zeros([n_noise,1]);
                            -Rtilde*dyn_model.gamma(dt), Rtilde, bhat_prev];

                        Ahouse = householder(A);
                        Rbar = Ahouse(n_noise+1:end, n_noise+1:end-1);
                        bbar = Ahouse(n_noise+1:end, end);
                    else
                        Rbar = Rhat_prev/Phi;
                        if triu_rbar
                            Rbar = householder(Rbar);
                        end
                        bbar = Rbar*dxbar;
                    end
                    
                end
                
        
                % Measurement Update
                % Identify visible stations
                stations = find(~isnan(Ydata(1,i,:)));
                if ~isempty(stations)
                    Vaug = kron(eye(length(stations)), chol(meas_model.R));
                    
                    % Pre-fit measurement residual
                    Y = reshape(Ydata(:,i,stations), [2*length(stations),1]);
                    Ystar = meas_model.measure(tdata(i), Xstar, stations);
                    dy = (Y - Ystar);
                    dYpre(:,i,stations) = reshape(dy, [2,1,length(stations)]);

                    % Whiten measurement partials and residuals
                    Htilde = meas_model.Htilde(tdata(i), Xstar, stations);
                    H_white = Vaug\Htilde;
                    dy_white = Vaug\dy;
                    
                    A = [Rbar, bbar;
                        H_white, dy_white];
                    A_house = householder(A);
                    Rhat = A_house(1:n_state,1:n_state);
                    bhat = A_house(1:n_state,end);
                    
                    inv_Rhat = inv(Rhat);
                    dxhat = inv_Rhat*bhat;

                    % Post-fit residuals
                    dYpost(:,i,stations) = reshape(dy - Htilde* dxhat, [2,1,length(stations)]);
              
                else
                    bhat = bbar;
                    dxhat = dxbar;
                    Rhat = Rbar;
                    inv_Rhat = inv(Rbar);
                end
                
                % Record estimates and update for next time step
                Xhist(:,i) = Xstar + dxhat;
                Phist(:,:,i) = inv_Rhat*inv_Rhat';
                bhat_prev = bhat;
                dxhat_prev = dxhat;
                Rhat_prev = Rhat;
            end
    
        end
    end

end