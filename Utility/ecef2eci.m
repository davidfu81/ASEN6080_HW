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