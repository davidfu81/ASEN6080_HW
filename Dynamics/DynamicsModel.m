classdef (Abstract) DynamicsModel
    % Base class for dynamics models
    
    properties
        n_state

        dyn_params
    end
    
    methods (Abstract)
        % Integrate equations of motion
        integrate_eom(obj, teval, X0)

        % Integrate equations of motion with STM
        integrate_eomwPhi(obj, teval, X0)

    end
end

