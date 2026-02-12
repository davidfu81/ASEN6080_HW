function A_matrix = dfdx_wJ2J3(r, v, mu, J2, J3)
% Outputs 9x9 for a 9d state [r,v, mu, J2, J3]
    RE = 6378;
    x = r(1); y = r(2); z = r(3);
    r = norm(r);

    A_matrix = zeros(9);

    % drdot/dv = I
    A_matrix(1:3,4:6) = eye(3);

    % dax/dx
    A_matrix(4,1) = mu*(3*x^2/r^5 - 1/r^3) + ...
        3/2*J2*RE^2*mu/r^5 * (5*z^2/r^2 - 1 + x^2*(5/r^2 - 35*z^2/r^4)) + ...
        1/2*mu*J3*RE^3/r^7 * (35*z^3/r^2 - 15*z + x^2*(105*z/r^2 - 315*z^3/r^4));
    % dax/dy = day/dx
    A_matrix(4,2) = 3*mu*x*y/r^5 + ...
        3/2*J2*RE^2*mu * (-35*z^2*y/r^9 + 5*y/r^7)*x + ...
        1/2*mu*J3*RE^3 * (-315*z^3*y/r^11 + 105*z*y/r^9)*x;
    A_matrix(5,1) = A_matrix(4,2);
    % dax/dz = daz/dx
    A_matrix(4,3) = 3*mu*x*z/r^5 + ...
        15/2*J2*RE^2*mu/r^7 * (3*z - 7*z^3/r^2)*x + ...
        15/2*mu*J3*RE^3/r^7 * (14*z^2/r^2 - 21*z^4/r^4 - 1)*x;
    A_matrix(6,1) = A_matrix(4,3);
    %day/dy
    A_matrix(5,2) = mu*(3*y^2/r^5 - 1/r^3) + ...
        3/2*J2*RE^2*mu/r^5 * (5*z^2/r^2 - 1 + y^2*(5/r^2 - 35*z^2/r^4)) + ...
        1/2*mu*J3*RE^3/r^7 * (35*z^3/r^2 - 15*z + y^2*(105*z/r^2 - 315*z^3/r^4));
    % day/dz = daz/dy
    A_matrix(5,3) = 3*mu*y*z/r^5 + ...
        15/2*J2*RE^2*mu/r^7 * (3*z - 7*z^3/r^2)*y + ...
        15/2*mu*J3*RE^3/r^7 * (14*z^2/r^2 - 21*z^4/r^4 - 1)*y;
    A_matrix(6,2) = A_matrix(5,3);
    % daz/dz
    A_matrix(6,3) = mu*(3*z^2/r^5 - 1/r^3) + ...
        3/2*J2*RE^2*mu/r^5 * (30*z^2/r^2 - 35*z^4/r^4 - 3) + ...
        5/2*mu*J3*RE^3/r^7 * (70*z^3/r^2 - 63*z^5/r^4 - 15*z);

    %da/dmu
    A_matrix(4,7) = -x/r^3 + 3/2*J2*RE^2*(5*z^2/r^7 - 1/r^5)*x + ...
        1/2*J3*RE^3*(35*z^3/r^9 - 15*z/r^7)*x;
    A_matrix(5,7) = -y/r^3 + 3/2*J2*RE^2*(5*z^2/r^7 - 1/r^5)*y + ...
        1/2*J3*RE^3*(35*z^3/r^9 - 15*z/r^7)*y;
    A_matrix(6,7) = -z/r^3 + 3/2*J2*RE^2*(5*z^2/r^7 - 3/r^5)*z + ...
        1/2*J3*RE^3*(35*z^4/r^9 - 30*z^2/r^7 + 3/r^5);

    %da/dJ2
    A_matrix(4,8) = 3/2*mu*RE^2*(5*z^2/r^7 - 1/r^5)*x;
    A_matrix(5,8) = 3/2*mu*RE^2*(5*z^2/r^7 - 1/r^5)*y;
    A_matrix(6,8) = 3/2*mu*RE^2*(5*z^2/r^7 - 3/r^5)*z;

    %da/dJ3
    A_matrix(4,9) = 1/2*mu*RE^3*(35*z^3/r^9 - 15*z/r^7)*x;
    A_matrix(5,9) = 1/2*mu*RE^3*(35*z^3/r^9 - 15*z/r^7)*y;
    A_matrix(6,9) = 1/2*mu*RE^3*(35*z^4/r^9 - 30*z^2/r^7 + 3/r^5);

end