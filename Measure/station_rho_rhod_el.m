function rho_rhod_el = station_rho_rhod_el(ref_time, ref_state, ...
    station_lat_deg, station_long_deg, theta0_deg)

    RE = 6378;
    wE = 360/(86400); %deg/s
    rho_rhod_el = nan([3, length(ref_time)]);

    % initial station ecef
    stationz = RE*sind(station_lat_deg);
    stationx = RE*cosd(station_lat_deg)*cosd(station_long_deg);
    stationy = RE*cosd(station_lat_deg)*sind(station_long_deg);
    station_ecef = [stationx; stationy; stationz];
    
    % Get measurements at each time
    for i = 1:length(ref_time)
        % Compute station ECI state
        theta = theta0_deg + wE*ref_time(i);
        station_eci = [cosd(theta), sind(theta), 0;
                        -sind(theta), cosd(theta), 0;
                        0, 0, 1]*station_ecef;
        station_vel = cross([0;0;deg2rad(wE)], station_eci);
        
        [rho, rhod, el] = G_rho_rhod_el(ref_state(:,i), station_eci, station_vel);

        if el >= 10
            rho_rhod_el(:,i) = [rho; rhod; el];
        end

    end
end

