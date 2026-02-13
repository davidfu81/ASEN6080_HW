classdef MeasurementModel

    properties 
        % Nonlinear measurement model
        G_func
        
        % Partials matrix for measurement model
        Htilde_func
        
        % Parameters for measurements
        meas_params

        % Number of scalars expected in state vector
        n_state

    end

    methods (Access = public)
        function obj = MeasurementModel(G_func, Htilde_func, meas_params, n_state)
            obj.G_func = G_func;
            obj.Htilde_func = Htilde_func;
            obj.meas_params = meas_params;
            obj.n_state = n_state;
        end
        function Xsim = integrate_eom(obj, teval, X0)
            if length(X0) ~= obj.n_state
                error("Input state must be length %d", obj.n_state)
            end

            opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);
            sol = ode45(@(t,y) obj.ode(t,y, obj.dyn_params), ...
                [teval(1), teval(end)], X0, opts);
            Xsim = deval(sol, teval);
        end

        function [Xsim, Phi] = integrate_eomwPhi(obj, teval, X0)
            if length(X0) ~= obj.n_state
                error("Input state must be length %d", obj.n_state)
            end
        
            opts = odeset('RelTol',1e-12, 'AbsTol', 1e-9);

            Phi0 = eye(obj.n_state);
            sol = ode45(@(t,y) obj.odewPhi(t,y), ...
                [teval(1), teval(end)], [X0; reshape(Phi0, [obj.n_state^2,1])], opts);
            Xaug = deval(sol, teval);
            Xsim = Xaug(1:obj.n_state,:);
        
            Phi = zeros(obj.n_state,obj.n_state, length(teval));
            for i = 1:length(teval)
                Phi(:,:,i) = reshape(Xaug(obj.n_state+1:end,i), [obj.n_state,obj.n_state]);
            end
        end
    end


end