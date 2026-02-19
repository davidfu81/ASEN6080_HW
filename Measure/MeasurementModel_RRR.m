classdef MeasurementModel_RRR < MeasurementModel

    properties 
        % ECEF positions of stations
        Recef_stations
    end

    methods (Access = public)
        function obj = MeasurementModel_RRR(Recef_stations, R)
            % Recef_stations = [Rs1, Rs2, ...] in ECEF
            obj.Recef_stations = Recef_stations;
            obj.R = R;
        end
        
        function [Y_rho_rhod, Y_el] = measure(obj, teval, Xeval, stations)
            %%%%
            % Nonlinear measurement model (for one state)
            % X = [x; y; z; vx; vy; vz; ...]
            % Rsi = [Rs1, Rs2, ...] in ECEF
            %%%%

            num_station = length(stations);
            Y_rho_rhod = zeros([2*num_station, 1]);
            Y_el = zeros([num_station, 1]);
            for k = 1:num_station

                [Rs, Vs] = CoordinateConverter.ecef2eci(teval, obj.Recef_stations(:,stations(k)));
            
                % Range and range rate
                rel_pos = Xeval(1:3) - Rs;
                rel_vel = Xeval(4:6) - Vs;
            
                rho = sqrt(sum(rel_pos.^2));
                rhod = sum(rel_pos.*rel_vel)./rho;
            
                % Elevation
                up_unit = Rs/Constants.RE;
                el = rad2deg(asin(sum(rel_pos.*up_unit)./rho));
           
                
                Y_rho_rhod(2*k-1:2*k,:) = [rho; rhod];
                Y_el(k) = el;
            end
            
        end

        function Htilde = Htilde(obj, teval, Xeval, stations)
            %%%%
            % X = [x; y; z; vx; vy; vz; ...]
            % Rsi = [Rs1, Rs2, ...] in ECEF
            %%%%
            num_station = length(stations);
            Htilde = zeros([2*num_station, length(Xeval)]);
            for k = 1:num_station
            
                [Rs, Vs] = CoordinateConverter.ecef2eci(teval, obj.Recef_stations(:,stations(k)));
            

                r = Xeval(1:3);
                v = Xeval(4:6);
                rho = norm(r-Rs);
                rhodot = (r-Rs)'*(v-Vs)/rho;
                
                % State partials
                Htilde(2*k-1, 1:3) = (r-Rs)'/rho;
                Htilde(2*k, 1:3) = (v-Vs)'/rho - rhodot/rho^2*(r-Rs)';
                Htilde(2*k, 4:6) = (r-Rs)'/rho;

           
            end
        end

        function Y = simulate_measure(obj, teval, Xeval, min_el)
            %%%%
            % Simulate nonlinear measurements for a whole trajectory
            % X = [x; y; z; vx; vy; vz; ... ]
            % R = measurement noise matrix 
            % min_el = minimum elevation (deg)
            %%%% 
            
            num_stations = size(obj.Recef_stations, 2);
            Y = nan(2, length(teval), num_stations);
            
            for i = 1:num_stations
                for k = 1:length(teval)
                    [Y_rho_rhod, Y_el] = obj.measure(teval(k), ...
                        Xeval(:,k), i);
                    if Y_el >= min_el
                        Y(:,k,i) = Y_rho_rhod + mvnrnd(zeros([1,2]), obj.R)';
                    end
                end
            end
        end
    end

end