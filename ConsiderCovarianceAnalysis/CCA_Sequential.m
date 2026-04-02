classdef CCA_Sequential

    methods
        function obj = CCA_Sequential()

        end

        function [Xhist, Phist, Xhist_c, Phist_c, Phist_xc, dYpre, dYpost, Psi_tot] = run_cca(~, tdata, Ydata, Xref0, Phat0, dyn_model, meas_model, consider_params)
            % Initialize
            n_state = length(Xref0);
            n_param = length(consider_params);

            Xhist = zeros([n_state,xqlength(tdata)]);
            Xhist_c = zeros([n_state,length(tdata)]);
            Phist = zeros([n_state,n_state,length(tdata)]);
            Phist_c = zeros([n_state,n_state,length(tdata)]);
            Phist_xc = zeros([n_state, n_param, length(tdata)]);
            dYpre = nan([size(Ydata)]);
            dYpost = nan([size(Ydata)]);
            
            X0 = Xref0;
            dx0 = zeros([n_state, 1]);
            dc = [consider_params.dc];
            Pcc = diag([consider_params.sigma^2]);
            param_names = [consider_params.param];
            
            
            % Iterate through observations
            for i = 1:length(tdata)
                % Time Update
                if i == 1
                    dxbar = dx0;
                    Pbar = Phat0;
                    Xstar = X0;

                    Htilde_x = meas_model.Htilde(0, X0, 1);
                    Htilde_c = meas_model.Htilde_c(0, X0, 1, param_names);
                    Sbar = -Pbar*Htilde_x'/(Htilde_x*Pbar*Htilde_x' + meas_model.R)*Htilde_c;
                    Pbar_c = Pbar;
                    Pbar_xc = Sbar*Pcc;
                    dxbar_c = dxbar;

                    Psi_tot = eye(n_state+n_param);
                else
                    % Integrate refrence trajectory with STM
                    [Xref, Phiref, Thetaref, Psi_ref] = dyn_model.integrate_eomwPhiTheta([tdata(i-1:i)], Xstar_prev, param_names);
                    
                    Phi = Phiref(:,:,end);
                    Theta = Thetaref(:,:,end);

                    Psi_tot = Psi_ref(:,:,end)*Psi_tot;
                    Xstar = Xref(:,end);
                    dt = tdata(i) - tdata(i-1);
                    dxbar = Phi*dxhat_prev;
                    Pbar = Phi*Phat_prev*Phi' + dyn_model.process_noise_covariance(dt, Xstar_prev);

                    % Consider time update
                    Sbar = Phi*S_prev + Theta;
                    dxbar_c = dxbar + Sbar*dc;
                    Pbar_c = Pbar + Sbar*Pcc*Sbar';
                    Pbar_xc = Sbar*Pcc;
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
                    dxhat = dxbar + K*(dy - Htilde*dxbar);
        
                    M = (eye(n_state) - K*Htilde);
                    Phat = M*Pbar*M' + K*Raug*K';
        
                    % Post-fit residuals
                    dYpost(:,i,stations) = reshape(dy - Htilde* dxhat, [2,1,length(stations)]);

                    % Consider measurement update
                    Htilde_c = meas_model.Htilde_c(tdata(i), Xstar, stations, param_names);
                    Shat = (eye(n_state) - K*Htilde)*Sbar - K*Htilde_c;
                    dxhat_c = dxhat + Shat * dc;
                    Phat_c = Phat + Shat*Pcc*Shat';
                    Phat_xc = Shat*Pcc;
                else
                    dxhat = dxbar;
                    Phat = Pbar;
                    dxhat_c = dxbar_c;
                    Phat_c = Pbar_c;
                    Shat = Sbar;
                    Phat_xc = Pbar_xc;
                end
                
                % Record estimates and update for next time step
                Xhist(:,i) = Xstar + dxhat;
                Xhist_c(:,i) = Xstar + dxhat_c;
                Phist(:,:,i) = Phat;
                Phist_c(:,:,i) = Phat_c;
                Phist_xc(:,:,i) = Phat_xc;

                dxhat_prev = dxhat;
                Phat_prev = Phat;
                Xstar_prev = Xstar;
                S_prev = Shat;
            end
    
        end
    end

end