classdef DynamicsModel_SNC < DynamicsModel
    % Dynamics model that adds SNC to another DynamicsModel
    % Note only diagonal Q matrices are supported
    properties
        base_dyn_model
        frame
        Q_pn
    end
    methods (Access = public)
        function obj = DynamicsModel_SNC(base_dyn_model, sig_x, sig_y, sig_z, frame)
            if nargin < 5
                frame = "ECI";
            end
            obj.base_dyn_model = base_dyn_model;
            obj.Q_pn = diag([sig_x^2, sig_y^2, sig_z^2]);
            obj.frame = frame;
        end
        function Xsim = integrate_eom(obj, teval, X0)
            Xsim = obj.base_dyn_model.integrate_eom(teval, X0);
        end

        function [Xsim, Phi] = integrate_eomwPhi(obj, teval, X0)
            [Xsim, Phi] = obj.base_dyn_model.integrate_eomwPhi(teval, X0);
        end

        function Q = process_noise_covariance(obj, dt, X)
            if obj.frame == "RIC"
                C_ric2eci = ric2eci(X);
                Q_pn_eci = C_ric2eci*obj.Q_pn*C_ric2eci';
            elseif obj.frame == "ECI"
                Q_pn_eci = obj.Q_pn;
            else
                error("Unrecognized SNC frame %s", obj.frame)
            end
            Q = dt^2*[dt^2/4*Q_pn_eci, dt/2*Q_pn_eci;
                dt/2*Q_pn_eci, Q_pn_eci];
        end

        function obj = set_Q(obj, sig_x, sig_y, sig_z, frame)
            if nargin < 5
                frame = obj.frame;
            end
            obj.Q_pn = diag([sig_x^2, sig_y^2, sig_z^2]);
            obj.frame = frame;
        end

        function dXdt = eom(obj, t, X)
            dXdt = obj.base_dyn_model.eom(t, X);
        end

        function Amatrix = dfdx(obj, X)
            Amatrix = obj.base_dyn_model.dfdx(X);
        end
    end


end