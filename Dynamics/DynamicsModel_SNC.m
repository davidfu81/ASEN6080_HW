classdef DynamicsModel_SNC < DynamicsModel
    % Dynamics model that adds SNC to another DynamicsModel
    % Note only diagonal Q matrices are supported
    properties
        base_dyn_model

        Q_pn
    end
    methods (Access = public)
        function obj = DynamicsModel_SNC(base_dyn_model, sig_x, sig_y, sig_z)
            obj.base_dyn_model = base_dyn_model;
            obj.Q_pn = diag([sig_x^2, sig_y^2, sig_z^2]);
        end
        function Xsim = integrate_eom(obj, teval, X0)
            Xsim = obj.base_dyn_model.integrate_eom(teval, X0);
        end

        function [Xsim, Phi] = integrate_eomwPhi(obj, teval, X0)
            [Xsim, Phi] = obj.base_dyn_model.integrate_eomwPhi(teval, X0);
        end

        function Q = process_noise_covariance(obj, dt)
            Q = dt^2*[dt^2/4*obj.Q_pn, dt/2*obj.Q_pn;
                dt/2*obj.Q_pn, obj.Q_pn];
        end

        function obj = set_Q(obj, sig_x, sig_y, sig_z)
            obj.Q_pn = diag([sig_x^2, sig_y^2, sig_z^2]);
        end
    end


end