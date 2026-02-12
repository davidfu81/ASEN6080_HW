% Utility function to plot measurement residuals
function [fig_resid, ax_resid] = plot_meas_resid(tdata, yest, ytrue, title)
    
    fig_resid = figure;
    ax_resid = [];
    set(gcf, "Position", [100 100 750 600])
    tl = tiledlayout(2,1);
    tl.TileSpacing = 'compact';
    tl.Padding = 'loose';
    
    ylabels = ["Rho Error [m]", "Rho Dot Error [m/s]"];
    for i = 1:2
        ax_resid(i) = nexttile;
        hold on; grid on;
        ax = gca;
        ax.XAxis.Exponent = 0;
        ax.YAxis.Exponent = 0;

        for j = 1:size(ytrue, 3)
            errors = ytrue(i,:,j) - yest(i,:,j);
            scatter(tdata, errors*1000, 4)
        end
        ylabel(ylabels(i))
    end
    xlabel("Time [s]")
    lgd_labels = strings([1, size(ytrue, 3)]);
    for i = 1:3
        lgd_labels(i) = sprintf("Station %d", i);
    end
    lgd = legend(lgd_labels, 'NumColumns', length(lgd_labels));
    lgd.Layout.Tile = 'north';
    sgtitle(title);
end