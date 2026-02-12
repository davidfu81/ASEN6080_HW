% Utility function to plot state estimation errors
function [fig_error, ax_error] = plot_filter_error(tdata, xest, Pest, xtrue, title)
    
    fig_error = figure;
    ax_error = [];
    set(gcf, "Position", [100 100 750 600])
    tl = tiledlayout(6,1);
    tl.TileSpacing = 'compact';
    tl.Padding = 'loose';
    
    ylabels = ["X Error [m]", "Y Error [m]", "Z Error [m]", ...
        "Vx Error [m/s]", "Vy Error [m/s]", "Vz Error [m/s]"];
    for i = 1:6
        ax_error(i) = nexttile;
        hold on; grid on;
        ax = gca;
        ax.XAxis.Exponent = 0;
        ax.YAxis.Exponent = 0;
        errorsi = xtrue(i,:) - xest(i,:);
        plot(tdata, errorsi*1000, 'LineStyle', '-', ...
            'Color', 'blue')
        sigma = reshape(sqrt(Pest(i,i,:)), 1, []);
        plot(tdata, -3*sigma*1000, 'LineStyle', '--', ...
            'Color', 'red');
        plot(tdata, zeros(size(tdata)), 'LineStyle', ':', ...
            'Color', 'black')
        plot(tdata, 3*sigma*1000, 'LineStyle', '--', ...
            'Color', 'red');
        ylabel(ylabels(i))
        % ylim([-max(abs(errorsi))*1.25, max(abs(errorsi))*1.25])
    end
    xlabel("Time [s]")
    lgd = legend(["Filter Error", "Filter 3\sigma Bounds", "Truth"], 'NumColumns', 3);
    lgd.Layout.Tile = 'north';
    sgtitle(title);
end