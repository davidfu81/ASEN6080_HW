classdef (Abstract) MeasurementModel

    properties

        % Measurement noise
        R

    end

    methods (Abstract)
        measure(obj, teval, Xeval, R)

        Htilde(obj, teval, Xeval)

        simulate_measure(obj, teval, Xeval, R)
    end


end