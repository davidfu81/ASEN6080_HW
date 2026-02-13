classdef (Abstract) MeasurementModel

    properties

        % Number of scalars expected in state vector
        n_state

    end

    methods (Abstract)
        measure(obj, teval, Xeval, R)

        Htilde(obj, teval, Xeval)

        simulate_measure(obj, teval, Xeval, R)
    end


end