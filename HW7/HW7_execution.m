clear
load("HW6/data_J2J3.mat", "tdata", "Ydata", "X0", "Xtrue", "dx0", "Rsi", ...
    "Phat0", "R", "period")

dyn_model_J2 = DynamicsModel_J2J3(struct('mu', Constants.MU, 'J2', ...
    Constants.J2, 'J3', 0, 'RE', Constants.RE), 6);
meas_model = MeasurementModel_RRR(Rsi, R);

cca = CCA_Sequential();
Phat0 = 1e4*eye(6)/1e6;
consider_params = struct('param', "J3", 'dc', Constants.J3, 'sigma', Constants.J3);
[Xhist, Phist, Xhist_c, Phist_c, Phist_xc, dYpre, dYpost] = cca.run_cca(tdata, Ydata, X0+dx0, Phat0, dyn_model_J2, meas_model, consider_params);

%%
close all
[fig_cca, ~, percent_bound] = plot_cca(tdata, Xhist*1e3, Phist*1e6, Phist_c*1e6, Xtrue*1e3, "J_3 Consider Covariance Analysis");
exportgraphics(fig_cca, "Figures/HW7/cca.png")
%%
[Xref, ~, ~, Psi_ref] = dyn_model_J2.integrate_eomwPhiTheta([tdata(1), tdata(end)], X0+dx0, "J3");

Psi_inv = inv(Psi_ref(:,:,end));
dxhatf = Xhist_c(:,end) - Xref(:,end);
dxhat0 = Psi_inv*[dxhatf; consider_params.dc];
Phatf = [Phist_c(:,:,end), Phist_xc(:,:,end); Phist_xc(:,:,end)', diag(consider_params.sigma^2)];
% Phatf = blkdiag(Phist_c(:,:,end), diag(consider_params.sigma^2));
Phat0 = Psi_inv*Phatf*Psi_inv';

dyn_model_J2J3 = DynamicsModel_J2J3(struct('mu', Constants.MU, 'J2', ...
    Constants.J2, 'J3', Constants.J3, 'RE', Constants.RE), 6);
[Xback, Phi_back, ~, Psi_back] = dyn_model_J2J3.integrate_eomwPhiTheta(tdata, Xref(:,1)+dxhat0(1:6), "J3");
Pprop = zeros(6,6,length(tdata));

for i = 1:length(tdata)
    Paug = Psi_back(:,:,i)*Phat0*Psi_back(:,:,i)';
    Pprop(:,:,i) = Paug(1:6,1:6);
    % Pprop(:,:,i) = Phi_back(:,:,i)*Phat0(1:6,1:6)*Phi_back(:,:,i)';
end

%%
fig_err = plot_cca(tdata, Xback*1e3, Pprop*1e6, [], Xtrue*1e3, "J_3 Consider State Errors Incorporating Full Dataset", false);
exportgraphics(fig_err, "Figures/HW7/err.png")
