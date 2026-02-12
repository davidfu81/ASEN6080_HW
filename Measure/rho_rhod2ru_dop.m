function ru_dop_el = rho_rhod2ru_dop(rho_rhod_el, ft)
    c = 2.998e5; % km/s

    ru_dop_el = nan(size(rho_rhod_el));

    ru_dop_el(1,:,:) = 221/749*rho_rhod_el(1,:,:)/c*ft;
    ru_dop_el(2,:,:) = -2*rho_rhod_el(2,:,:)/c*ft;
    ru_dop_el(3,:,:) = rho_rhod_el(3,:,:);
end

