classdef Constants

    properties (Constant)
        % Earth properties
        RE = 6378136.3/1e3; %km
        MU = 3.986004415e14/1e9; %km3/s2
        J2 = 1.082626925638815e-3;
        WE = 7.2921158553e-5; % rad/s
        WE_VEC = [0;0;7.2921158553e-5];

        % Atmosphere Model
        RHO0 = 3.614e-13*1e9 %kg/km3
        R0 = (700e3+6378136.3)/1e3; % km
        H = 88667/1e3 % km

    end

    methods (Access = private)
        function obj = Constants()
        end
    end


end