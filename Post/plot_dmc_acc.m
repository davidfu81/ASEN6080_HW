function fig_dmc_acc = plot_dmc_acc(tdata, acc_est, acc_true, title_str)

    ylabels = ["Ax Error [mm/s^2]", "Ay Error [mm/s^2]", "Ax Error [mm/s^2]"];
    fig_dmc_acc = figure;
    set(gcf, "Position", [100 100 750 500])
    tl = tiledlayout(3,1); tl.TileSpacing = 'compact'; tl.Padding = 'loose';
    for i = 1:3
        nexttile; grid on; hold on;
        plot(tdata, acc_true(i,:)*1e6)
        plot(tdata, acc_est(i,:)*1e6)
        ylabel(ylabels(i))
        ax = gca;
        ax.XAxis.Exponent = 0;
        ax.YAxis.Exponent = 0;
    end
    xlabel("Time [s]")
    lgd = legend(["J3 Acceleration", "DMC Acceleration Estimate"], 'NumColumns', 2);
    lgd.Layout.Tile = 'north';
    sgtitle(title_str)
end