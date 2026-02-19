% Utility function to plot state estimation errors
function [fig_error, ax_error] = plot_filter_error(tdata, xest, Pest, xtrue, title, ylabels)
    if nargin < 6
        ylabels = ["X Error [m]", "Y Error [m]", "Z Error [m]", ...
        "Vx Error [m/s]", "Vy Error [m/s]", "Vz Error [m/s]"];
    end
    fig_error = figure;
    ax_error = [];
    set(gcf, "Position", [100 100 750 400+100*size(xest,1)/3])
    tl = tiledlayout(size(xest,1),1);
    tl.TileSpacing = 'compact';
    tl.Padding = 'loose';
    
    for i = 1:size(xest,1)
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
        if max(3*sigma) > max(abs(errorsi))*50
            ylim([-max(abs(errorsi)*1e3)*1.25, max(abs(errorsi)*1e3)*1.25])
        end
    end
    xlabel("Time [s]")
    lgd = legend(["Filter Error", "Filter 3\sigma Bounds", "Truth"], 'NumColumns', 3);
    lgd.Layout.Tile = 'north';
    sgtitle(title);
end