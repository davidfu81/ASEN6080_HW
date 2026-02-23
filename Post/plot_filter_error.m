% Utility function to plot state estimation errors
function [fig_error, ax_error] = plot_filter_error(tdata, xest, Pest, xtrue, title_str, rms_mask, components, units)
   if nargin < 6
        rms_mask = true(size(tdata));
    end
    if nargin < 7
        components = ["X", "Y", "Z", "Vx", "Vy", "Vz"];
        units = ["m", "m", "m", "m/s", "m/s", "m/s"];
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
        plot(tdata, errorsi, 'LineStyle', '-', ...
            'Color', 'blue')
        rms_error = sqrt(mean(errorsi(rms_mask).^2));
        title(sprintf("RMS Error: %.3e %s", rms_error, units(i)))

        sigma = reshape(sqrt(Pest(i,i,:)), 1, []);
        plot(tdata, -3*sigma, 'LineStyle', '--', ...
            'Color', 'red');
        plot(tdata, zeros(size(tdata)), 'LineStyle', ':', ...
            'Color', 'black')
        plot(tdata, 3*sigma, 'LineStyle', '--', ...
            'Color', 'red');
        ylabel(sprintf("%s Error [%s]", components(i), units(i)))
        if max(3*sigma) > max(abs(errorsi))*50
            ylim([-max(abs(errorsi))*1.25, max(abs(errorsi))*1.25])
        end
    end
    xlabel("Time [s]")
    lgd = legend(["Filter Error", "Filter 3\sigma Bounds", "Truth"], 'NumColumns', 3);
    lgd.Layout.Tile = 'north';
    sgtitle(title_str);
end