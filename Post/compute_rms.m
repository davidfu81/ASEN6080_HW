function [rms_state_comp, rms_state_3d, rms_meas] = compute_rms(Xtrue, Xest, Ytrue, Yest)
    
    dX = Xtrue - Xest;
    dY = Ytrue - Yest;
    dY = dY(:, ~all(isnan(dY(1,:,:)),3),:);
    dY(isnan(dY)) = 0;

    rms_state_comp = sqrt(1/size(Xtrue, 2)*sum(dX.^2, 2));

    rms_state_3d(1) = sqrt(1/size(Xtrue, 2)*sum(sqrt(sum(dX(1:3,:).^2)).^2));
    rms_state_3d(2) = sqrt(1/size(Xtrue, 2)*sum(sqrt(sum(dX(4:6,:).^2)).^2));

    rms_meas = sqrt(1/size(dY, 2)*sum(sum(dY.^2, 2),3));
end