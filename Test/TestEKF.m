classdef TestEKF < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function testEKF(testCase)
            close all
            load("Test/test_ekf.mat",  "dYpost", "dYpre", "Pest", "Phat0", ...
                "R", "Rsi", "Xest", "Ydata", "X0", "teval");
            addpath("./..")

            
            params = struct('mu', Constants.MU, 'J2', Constants.J2, 'J3', 0, 'RE', Constants.RE);
            
            dyn_model = DynamicsModel_J2J3(params, 6);
            meas_model = MeasurementModel_RRR(Rsi, R);

            ekf = Filter_EKF_slow();
            
            [Xest_out, Pest_out, dYpre_out, dYpost_out] = ekf.run_filter(teval, Ydata, X0, Phat0, dyn_model, meas_model, 0);
            

            testCase.verifyEqual(Xest_out(1:3,:), Xest(1:3,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(Xest_out(4:6,:), Xest(4:6,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)

            testCase.verifyEqual(dYpre_out(1,:), dYpre(1,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(dYpre_out(2,:), dYpre(2,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)

            testCase.verifyEqual(dYpost_out(1,:), dYpost(1,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(dYpost_out(2,:), dYpost(2,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)
            
            % This tolerance was loosened because of numerical differences
            % depending on how long integrations are run for in the filter
            testCase.verifyEqual(Pest_out, Pest, 'RelTol', 3e-2)

        end

    end

end