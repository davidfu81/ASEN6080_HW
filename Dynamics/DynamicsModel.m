classdef (Abstract) DynamicsModel < handle
    % Base class for dynamics models
    
    
    methods (Abstract)
        % Integrate equations of motion
        integrate_eom(obj, teval, X0)

        % Integrate equations of motion with STM
        integrate_eomwPhi(obj, teval, X0)

        % Get process noise covariance inflation
        process_noise_covariance(obj) 

    end
end

