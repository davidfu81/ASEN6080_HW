classdef DynamicsModel_J2J3 < DynamicsModel
    % Dynamics model that includes J2 and J3 perturbations
    properties
        n_state

        dyn_params
        
    end
    methods (Access = public)
        function obj = DynamicsModel_J2J3(dyn_params, n_state)
            obj.n_state = n_state;
            obj.dyn_params = dyn_params;
        end
        function Xsim = integrate_eom(obj, teval, X0)
            if length(X0) ~= obj.n_state
                error("Input state must be length %d", obj.n_state)
            end

            opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);
            sol = ode45(@(t,y) obj.eom(t,y), [teval(1), teval(end)], X0, opts);
            Xsim = deval(sol, teval);
        end

        function [Xsim, Phi] = integrate_eomwPhi(obj, teval, X0)
            if length(X0) ~= obj.n_state
                error("Input state must be length %d", obj.n_state)
            end
        
            opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);

            Phi0 = eye(obj.n_state);
            sol = ode45(@(t,y) obj.eomwPhi(t,y), ...
                [teval(1), teval(end)], [X0; reshape(Phi0, [obj.n_state^2,1])], opts);
            Xaug = deval(sol, teval);
            Xsim = Xaug(1:obj.n_state,:);
        
            Phi = zeros(obj.n_state,obj.n_state, length(teval));
            for i = 1:length(teval)
                Phi(:,:,i) = reshape(Xaug(obj.n_state+1:end,i), [obj.n_state,obj.n_state]);
            end
        end

        function [Xsim, Phi, Theta, Psi] = integrate_eomwPhiTheta(obj, teval, X0, consider_params)
            n = obj.n_state;
            m = length(consider_params);

            opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);

            Phi0 = eye(n);
            Theta0 = zeros([n,1]);
            Xaug0 = [X0; reshape(Phi0, [n^2,1]); reshape(Theta0, [n*m,1])];
            sol = ode45(@(t,y) obj.eomwPhiTheta(t,y, consider_params), ...
                [teval(1), teval(end)], Xaug0, opts);
            Xaug = deval(sol, teval);
            Xsim = Xaug(1:n,:);
        
            Phi = zeros(n,n, length(teval));
            Theta = zeros(n,m, length(teval));
            Psi = zeros(n+m, n+m, length(teval));
            for i = 1:length(teval)
                Phi(:,:,i) = reshape(Xaug(n+1:n+n^2,i), [n,n]);
                Theta(:,:,i) = reshape(Xaug(n+n^2+1:n+n^2+n*m,i), [n,m]);
                Psi(:,:,i) = [Phi(:,:,i), Theta(:,:,i); zeros(m, n), eye(m)];
            end
        end

        function Q = process_noise_covariance(obj, ~, ~)
            Q = zeros(obj.n_state);
        end

        function dXdt = eom(obj, t, X)
            % X = [r; v]
            RE = obj.dyn_params.RE;
            mu = obj.dyn_params.mu;
            J2 = obj.dyn_params.J2;
            J3 = obj.dyn_params.J3;
            
            dXdt = zeros(size(X));
        
            % dr/dt = v
            dXdt(1:3) = X(4:6);
        
            % dv/dt = acceleration
            r = norm(X(1:3));
            z  = X(3);
            a2b = -mu/r^3 * X(1:3);
            aJ2 = 3/2*J2*RE^2*mu/r^5* ...
                [(5*z^2/r^2-1); (5*z^2/r^2-1); (5*z^2/r^2-3)].*X(1:3);
            aJ3 = 1/2*mu*J3*RE^3 * ...
                [(35*z^3/r^9-15*z/r^7)*X(1); ...
                (35*z^3/r^9-15*z/r^7)*X(2); ...
                (35*z^4/r^9-30*z^2/r^7+3/r^5)];
            dXdt(4:6) = a2b + aJ2 + aJ3;
            
        end

        function A_matrix = dfdx(obj, X)
        % Outputs 6x6 for a 6d state [x; y; z; vx; vy; vz]
            RE = obj.dyn_params.RE;
            mu = obj.dyn_params.mu;
            J2 = obj.dyn_params.J2;
            J3 = obj.dyn_params.J3;
            x = X(1); y = X(2); z = X(3);
            r = norm(X(1:3));
        
            A_matrix = zeros(length(X));
        
            % drdot/dv = I
            A_matrix(1:3,4:6) = eye(3);
        
            % dax/dx
            A_matrix(4,1) = mu*(3*x^2/r^5 - 1/r^3) + ...
                3/2*J2*RE^2*mu/r^5 * (5*z^2/r^2 - 1 + x^2*(5/r^2 - 35*z^2/r^4)) + ...
                1/2*mu*J3*RE^3/r^7 * (35*z^3/r^2 - 15*z + x^2*(105*z/r^2 - 315*z^3/r^4));
            % dax/dy = day/dx
            A_matrix(4,2) = 3*mu*x*y/r^5 + ...
                3/2*J2*RE^2*mu * (-35*z^2*y/r^9 + 5*y/r^7)*x + ...
                1/2*mu*J3*RE^3 * (-315*z^3*y/r^11 + 105*z*y/r^9)*x;
            A_matrix(5,1) = A_matrix(4,2);
            % dax/dz = daz/dx
            A_matrix(4,3) = 3*mu*x*z/r^5 + ...
                15/2*J2*RE^2*mu/r^7 * (3*z - 7*z^3/r^2)*x + ...
                15/2*mu*J3*RE^3/r^7 * (14*z^2/r^2 - 21*z^4/r^4 - 1)*x;
            A_matrix(6,1) = A_matrix(4,3);
            %day/dy
            A_matrix(5,2) = mu*(3*y^2/r^5 - 1/r^3) + ...
                3/2*J2*RE^2*mu/r^5 * (5*z^2/r^2 - 1 + y^2*(5/r^2 - 35*z^2/r^4)) + ...
                1/2*mu*J3*RE^3/r^7 * (35*z^3/r^2 - 15*z + y^2*(105*z/r^2 - 315*z^3/r^4));
            % day/dz = daz/dy
            A_matrix(5,3) = 3*mu*y*z/r^5 + ...
                15/2*J2*RE^2*mu/r^7 * (3*z - 7*z^3/r^2)*y + ...
                15/2*mu*J3*RE^3/r^7 * (14*z^2/r^2 - 21*z^4/r^4 - 1)*y;
            A_matrix(6,2) = A_matrix(5,3);
            % daz/dz
            A_matrix(6,3) = mu*(3*z^2/r^5 - 1/r^3) + ...
                3/2*J2*RE^2*mu/r^5 * (30*z^2/r^2 - 35*z^4/r^4 - 3) + ...
                5/2*mu*J3*RE^3/r^7 * (70*z^3/r^2 - 63*z^5/r^4 - 15*z);
        
        end

        function B_matrix = dfdc(obj, X, consider_params)
            
            RE = obj.dyn_params.RE;
            mu = obj.dyn_params.mu;
            z  = X(3);
            r = norm(X(1:3));

            B_matrix = zeros([length(X), length(consider_params)]);
            
            for i = 1:length(consider_params)
                if consider_params(i) == "J3"
                    %da/dJ3
                    B_matrix(4:6,i) = 1/2*mu*RE^3 * ...
                        [(35*z^3/r^9-15*z/r^7)*X(1); ...
                        (35*z^3/r^9-15*z/r^7)*X(2); ...
                        (35*z^4/r^9-30*z^2/r^7+3/r^5)];
                end
            end
            
        end

        function dXaugdt = eomwPhi(obj, t, Xaug)
            %%%%%%
            % Xaug is length (n^2 + n) where n is number of states
            %%%%%%
        
            dXaugdt = zeros(size(Xaug));
        
            dXaugdt(1:6) = obj.eom(t, Xaug(1:obj.n_state));

            phi = reshape(Xaug(obj.n_state+1:end), [obj.n_state, obj.n_state]);
            A = obj.dfdx(Xaug(1:obj.n_state));
            phidot = A*phi;
        
            dXaugdt(obj.n_state+1:end) = reshape(phidot, [obj.n_state^2,1]);
        end

        function dXaugdt = eomwPhiTheta(obj, t, Xaug, consider_params)
            n = obj.n_state;
            m = length(consider_params);
            dXaugdt = zeros(size(Xaug));
            
            % State evolution
            dXaugdt(1:n) = obj.eom(t, Xaug(1:n));
            
            % Phi evolution
            phi = reshape(Xaug(n+1:n^2+n), [n, n]);
            A = obj.dfdx(Xaug(1:n));
            phidot = A*phi;
        
            dXaugdt(n+1:n^2+n) = reshape(phidot, [n^2,1]);

            % Theta evolution
            B = obj.dfdc(Xaug(1:n), consider_params);
            theta = Xaug(n^2+n+1:n*m+n^2+n);
            thetadot = A*theta + B;
            dXaugdt(n^2+n+1:n*m+n^2+n) = reshape(thetadot, [n*m,1]);
        end

        
    end


end