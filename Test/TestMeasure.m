classdef TestMeasure < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function testMeasure(testCase)
            load("Test/test_measure.mat",  "X0", "dx0", "alpha", "Ylin", "Ynonlin", "Rsi");
            addpath("./..")

            
            measure_model = MeasurementModel_RRR(Rsi, 6);
            Ynonlin_out = zeros([2, length(alpha),3]);
            Ylin_out = Ynonlin_out;
            
            for i = 1:length(alpha)
                X = X0+alpha(i)*dx0;
                Yref = measure_model.measure(0, X0, [1,2,3], zeros(2));
                Htilde = measure_model.Htilde(0, X0, [1,2,3]);
                Ylin_out(:,i,:) = reshape(Yref + Htilde*alpha(i)*dx0, [2,1,3]);
            
                Ynonlin_out(:,i,:) = reshape(measure_model.measure(0, X, [1,2,3], zeros(2)), [2,1,3]);
            end

            testCase.verifyEqual(Ynonlin_out(1,:), Ynonlin(1,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(Ynonlin_out(2,:), Ynonlin(2,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)

            testCase.verifyEqual(Ylin_out(1,:), Ylin(1,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(Ylin_out(2,:), Ylin(2,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)

        end

    end

end