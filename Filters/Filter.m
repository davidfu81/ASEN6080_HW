classdef (Abstract) Filter

    properties
        % Dynamics Model Object
        dyn_model

        % Measurement Model Object
        meas_model
        
    end

    methods (Abstract)
        run_filter(tdata, Ydata, Xref0, Phat0)
        
    end

end