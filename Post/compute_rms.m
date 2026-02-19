function [rms_state_comp, rms_state_3d, rms_meas] = compute_rms(Xtrue, Xest, dYest)
    
    dX = Xtrue - Xest;
    dX = dX(:,~all(isnan(dYest(1,:,:)),3));

    rms_state_comp = sqrt(1/size(dX, 2)*sum(dX.^2, 2));

    rms_state_3d(1) = sqrt(1/size(dX, 2)*sum(sqrt(sum(dX(1:3,:).^2)).^2));
    rms_state_3d(2) = sqrt(1/size(dX, 2)*sum(sqrt(sum(dX(4:6,:).^2)).^2));

    rms_meas = sqrt(mean(dYest.^2, [2, 3], 'omitnan'));
end