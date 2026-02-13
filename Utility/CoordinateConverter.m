classdef CoordinateConverter
    % Utility class with static function to make conversions

    methods (Static)
        function [Reci, Veci] = ecef2eci(t, Recef)
            % Transform an ECEF position (assumed to be fixed) to its ECI position
            % and velocity
        
            theta = Constants.WE*t;
            C = [cos(theta), -sin(theta), 0;
                sin(theta), cos(theta), 0;
                0, 0, 1];
            Reci = C*Recef;
            Veci = cross(Constants.WE_VEC, Reci);
        end

        function C_ric2eci = C_ric2eci(Xeci)
            % Get the DCM to convert from ECI to RIC frame
            % X = [r; v]
        
            r_hat = Xeci(1:3)/norm(Xeci(1:3));
            c_hat = cross(Xeci(1:3), Xeci(4:6));
            c_hat = c_hat/norm(c_hat);
            i_hat = cross(c_hat, r_hat);
        
            C_ric2eci = [r_hat, i_hat, c_hat];
        end

        function cart = oscelt2cart(sma, ecc, inc, raan, argp, tanom, mu)
            % Convert from keplerian elements to cartesian coordinates
            % inc, raan, argp, tanom in degrees
        
            if ecc < 1e-8
	            e_hat = [cos(raan), sin(raan), 0];
	            h_hat = [sin(inc)*sin(raan), -sin(inc)*cos(raan), cos(inc)];
	            e_perp = cross(h_hat, e_hat);
            else
	            e_hat = [cos(argp)*cos(raan)-cos(inc)*sin(argp)*sin(raan),...
					            (cos(argp)*sin(raan) + cos(inc)*sin(argp)*cos(raan)),...
					            sin(argp)*sin(inc)];
	            e_perp = [-sin(argp)*cos(raan)-cos(inc)*cos(argp)*sin(raan),...
					            (-sin(argp)*sin(raan) + cos(inc)*cos(argp)*cos(raan)),...
					            cos(argp)*sin(inc)];
            end
        
            p = sma*(1-ecc^2);
            r = p/(1+ecc*cos(tanom))*(cos(tanom)*e_hat+sin(tanom)*e_perp);
            v = sqrt(mu/(sma*(1-ecc^2)))*(-sin(tanom)*e_hat+(ecc+cos(tanom))*e_perp);
        
            cart =  [r'; v'];
        end
    end


end

