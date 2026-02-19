classdef TestDynamics < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
        
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function testDynamicsJ2(testCase)
            load("Test/test_dynamics.mat",  "X0", "dx0", "teval", "Phi", "Xhist", "XhistPhi");
            addpath("./..")

            params = struct('mu', Constants.MU, 'J2', Constants.J2, 'J3', 0, 'RE', Constants.RE);
            
            model = DynamicsModel_J2J3(params, 6); 
            Xhist_out = model.integrate_eom(teval, X0+dx0);

            [Xref_out, Phi_out] = model.integrate_eomwPhi(teval, X0);
            
            XhistPhi_out = zeros([length(X0),length(teval)]);
            for i = 1:length(teval)
                XhistPhi_out(:,i) = Xref_out(:,i) + Phi_out(:,:,i)*dx0;
            end

            
            testCase.verifyEqual(Xhist_out(1:3,:), Xhist(1:3,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(Xhist_out(4:6,:), Xhist(4:6,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)
            
            testCase.verifyEqual(Phi_out, Phi, 'RelTol', 1e-9)
            testCase.verifyEqual(XhistPhi_out(1:3,:), XhistPhi(1:3,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(XhistPhi_out(4:6,:), XhistPhi(4:6,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)

        end

        function testDynamicsJ3(testCase)
            load("Test/test_dynamicsJ3.mat",  "X0", "dx0", "teval", "Phi", "Xhist", "XhistPhi");
            addpath("./..")

            params = struct('mu', Constants.MU, 'J2', Constants.J2, 'J3', Constants.J3, 'RE', Constants.RE);
            
            model = DynamicsModel_J2J3(params, 6); 
            Xhist_out = model.integrate_eom(teval, X0+dx0);

            [Xref_out, Phi_out] = model.integrate_eomwPhi(teval, X0);
            
            XhistPhi_out = zeros([length(X0),length(teval)]);
            for i = 1:length(teval)
                XhistPhi_out(:,i) = Xref_out(:,i) + Phi_out(:,:,i)*dx0;
            end

            
            testCase.verifyEqual(Xhist_out(1:3,:), Xhist(1:3,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(Xhist_out(4:6,:), Xhist(4:6,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)
            
            testCase.verifyEqual(Phi_out, Phi, 'RelTol', 1e-9)
            testCase.verifyEqual(XhistPhi_out(1:3,:), XhistPhi(1:3,:), 'RelTol', 1e-9, 'AbsTol', 1e-6)
            testCase.verifyEqual(XhistPhi_out(4:6,:), XhistPhi(4:6,:), 'RelTol', 1e-9, 'AbsTol', 1e-9)

        end

        function testHW1(testCase)
            addpath("./..")
            sol_1c = jsondecode(fileread('HW1/prob1c_solution.json'));
            state = sol_1c.inputs.state;
            params = struct('mu', state.mu, 'J2',state.J2, 'J3', state.J3, 'RE', 6378);
            
            model = DynamicsModel_J2J3(params, 6); 
            A_matrix = model.dfdx([state.r; state.v]);
            testCase.verifyEqual(A_matrix, sol_1c.outputs.A_matrix.values(1:6,1:6), 'RelTol', 2e-14)
        end
    end

end