classdef DynamicsModel_DMC < DynamicsModel
    % Dynamics model that adds DMC to another DynamicsModel
    properties
        base_model
        Q_dmc
        tau
    end
    methods (Access = public)
        function obj = DynamicsModel_DMC(base_model, sig_x, sig_y, sig_z, tau)

            obj.base_model = base_model;
            obj.Q_dmc = diag([sig_x^2, sig_y^2, sig_z^2]);
            obj.tau = tau;
        end
        function Xsim = integrate_eom(obj, teval, X0)
            opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);
            sol = ode45(@(t,y) obj.eom(t,y), [teval(1), teval(end)], X0, opts);
            Xsim = deval(sol, teval);
        end

        function [Xsim, Phi] = integrate_eomwPhi_alt(obj, teval, X0)
            opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);

            Phi0 = eye(9);
            sol = ode45(@(t,y) obj.eomwPhi_alt(t,y), ...
                [teval(1), teval(end)], [X0; reshape(Phi0, [9^2,1])], opts);
            Xaug = deval(sol, teval);
            Xsim = Xaug(1:9,:);
        
            Phi = zeros(9,9, length(teval));
            for i = 1:length(teval)
                Phi(:,:,i) = ...
                    reshape(Xaug(9+1:end,i), size(Phi0));
              
            end
        end

        function [Xsim, Phi] = integrate_eomwPhi(obj, teval, X0)
            opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);

            Phi0 = eye(6);
            sol = ode45(@(t,y) obj.eomwPhi(t,y), ...
                [teval(1), teval(end)], [X0; reshape(Phi0, [6^2,1])], opts);
            Xaug = deval(sol, teval);
            Xsim = Xaug(1:9,:);

            Phi = zeros(9,9, length(teval));
            for i = 1:length(teval)
                Phi(1:6,1:6,i) = ...
                    reshape(Xaug(9+1:end,i), size(Phi0));
                Phi(6+1:end, 1:6, i) = zeros([3, 6]);
                dt = (teval(i)-teval(1));
                exp_dt_tau = exp(-dt/obj.tau);
                Phi(1:3, 7:9,i) = (dt*obj.tau+obj.tau^2*(exp_dt_tau-1)) * eye(3);
                Phi(4:6, 7:9,i) = obj.tau*(1-exp_dt_tau)*eye(3);
                Phi(7:9, 7:9,i) = exp_dt_tau*eye(3);
            end
        end

        function Q = process_noise_covariance(obj, dt, X)
            Q = zeros(9);
            exp_dt_tau = exp(-(dt)/obj.tau);
            % Q_w_ii
            Q(1:3,1:3) = obj.Q_dmc*obj.tau^2 * (dt^3/3 - dt^2*obj.tau +...
                dt*obj.tau^2*(1 - 2*exp_dt_tau) + 1/2*obj.tau^3*(1-exp_dt_tau^2));
            %Q_w_ivi
            Q(1:3,4:6) = obj.Q_dmc*obj.tau^2 * (dt^2/2 - dt*obj.tau*(1-exp_dt_tau) + ...
                obj.tau^2*(1/2 - exp_dt_tau + 1/2*exp_dt_tau^2));
            Q(4:6,1:3) = Q(1:3,4:6);
            %Q_w_vivi
            Q(4:6, 4:6) = obj.Q_dmc*obj.tau^2 * (dt - obj.tau*(3/2 + 1/2*exp_dt_tau^2 - 2*exp_dt_tau));
            %Q_w_iwi
            Q(1:3,7:9) =  obj.Q_dmc*obj.tau^2 * (obj.tau/2*(1-exp_dt_tau^2) - dt*exp_dt_tau);
            Q(7:9,1:3) = Q(1:3,7:9);
            %Q_w_viwi
            Q(4:6,7:9) = obj.Q_dmc*obj.tau^2 * (1/2*(1+exp_dt_tau^2) - exp_dt_tau);
            Q(7:9,4:6) = Q(4:6, 7:9);
            %Q_w_wiwi
            Q(7:9,7:9) = obj.Q_dmc*obj.tau/2*(1-exp_dt_tau^2);
            
        end

        function obj = set_Q_dmc(obj, sig_x, sig_y, sig_z)
            obj.Q_dmc = diag([sig_x^2, sig_y^2, sig_z^2]);
        end

        function obj = set_tau(obj, tau)
            obj.tau = tau;
        end

        function dXdt = eom(obj, t, X)
            dXdt = zeros(size(X));
            dXdt(1:6) = obj.base_model.eom(t, X(1:6)) +...
                [zeros([3,1]); X(6+1:end)];
            dXdt(6+1:9) = -1/obj.tau*X(7:9);
        end
         
        function dXaugdt = eomwPhi(obj, t, Xaug)
            % Xaug is 9 + 6x6 states
            dXaugdt = zeros(size(Xaug));
            dXaugdt(1:9) = obj.eom(t, Xaug(1:9));
            Phi = reshape(Xaug(9+1:end), [6,6]);
            Phidot = obj.base_model.dfdx(Xaug(1:6))*Phi;
            dXaugdt(9+1:end) = reshape(Phidot, [6^2, 1]);
        end

        function dXaugdt = eomwPhi_alt(obj, t, Xaug)
            % Xaug is 9 + 9x9 states
            dXaugdt = zeros(size(Xaug));
            dXaugdt(1:9) = obj.eom(t, Xaug(1:9));
            Phi = reshape(Xaug(9+1:end), [9,9]);
            Phidot = obj.dfdx(Xaug(1:9))*Phi;
            dXaugdt(9+1:end) = reshape(Phidot, [9^2, 1]);
        end
        
        % Implemented because required but don't use for integration
        function Amatrix = dfdx(obj, X)
            Amatrix = zeros(9);
            Amatrix(1:6,1:6) = obj.base_model.dfdx(X(1:6));
            Amatrix(4:6,7:9) = eye(3);
            Amatrix(7:9,7:9) = -1/obj.tau*eye(3);
        end
    end


end