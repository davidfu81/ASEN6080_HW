% Utility function to plot measurement residuals
function [fig_resid, ax_resid] = plot_meas_resid(tdata, yresid, title_str, ...
    components, units, stations)
    if nargin < 4
        components = ["Range", "Range-Rate"];
        units = ["m", "m/s"];
    end
    if nargin < 6
        stations = [1, 2, 3];
    end
    fig_resid = figure;
    ax_resid = [];
    set(gcf, "Position", [100 100 750 500])
    tl = tiledlayout(length(components),1);
    tl.TileSpacing = 'compact';
    tl.Padding = 'loose';
    
    for i = 1:length(components)
        ax_resid(i) = nexttile;
        hold on; grid on;
        ax = gca;
        ax.XAxis.Exponent = 0;
        ax.YAxis.Exponent = 0;

        for j = 1:size(yresid, 3)
            scatter(tdata, yresid(i,:,j), 4)
        end
        ylabel(sprintf("%s [%s]", components(i), units(i)))
        rms_resid = sqrt(sum(sum(yresid(i,~isnan(yresid(i,:,:))).^2))/length(tdata));
        title(sprintf("RMS: %.3e %s", rms_resid, units(i) ))
    end
    xlabel("Time [s]")
    lgd_labels = strings([1, size(yresid, 3)]);
    for i = 1:3
        lgd_labels(i) = sprintf("Station %d", stations(i));
    end
    lgd = legend(lgd_labels, 'NumColumns', length(lgd_labels));
    lgd.Layout.Tile = 'north';
    sgtitle(title_str);
end