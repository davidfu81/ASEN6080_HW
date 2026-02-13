function C_ric2eci = ric2eci(Xeci)
    % X = [r; v]

    r_hat = Xeci(1:3)/norm(Xeci(1:3));
    c_hat = cross(Xeci(1:3), Xeci(4:6));
    c_hat = c_hat/norm(c_hat);
    i_hat = cross(c_hat, r_hat);

    C_ric2eci = [r_hat, i_hat, c_hat];
end