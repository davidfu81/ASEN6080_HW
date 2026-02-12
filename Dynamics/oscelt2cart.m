function cart = oscelt2cart(sma, ecc, inc, raan, argp, tanom, mu)
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