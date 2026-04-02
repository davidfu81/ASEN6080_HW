
function [fig_cca, ax_cca, percent_bound] = plot_cca(tdata, xest, Pest, Pc, xtrue, title_str, scale2err)
   % if nargin < 7
   %      rms_mask = true(size(tdata));
   %  end
   %  if nargin < 8
   %      components = ["X", "Y", "Z", "Vx", "Vy", "Vz"];
   %      units = ["m", "m", "m", "m/s", "m/s", "m/s"];
   %  end

    if nargin < 7
       scale2err = true;
    end
    rms_mask = true(size(tdata));
    components = ["X", "Y", "Z", "Vx", "Vy", "Vz"];
    units = ["m", "m", "m", "m/s", "m/s", "m/s"];
    percent_bound = zeros(6,2);
    
    fig_cca = figure;
    ax_cca = [];
    set(gcf, "Position", [100 100 750 400+100*size(xest,1)/3])
    tl = tiledlayout(size(xest,1),1);
    tl.TileSpacing = 'compact';
    tl.Padding = 'loose';
    
    for i = 1:size(xest,1)
        ax_cca(i) = nexttile;
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
    
        plot(tdata, -2*sigma, 'LineStyle', '--', ...
            'Color', 'red');
        percent_bound(i,1) = sum(abs(errorsi) < 2*sigma)/length(tdata);
        if ~isempty(Pc)
            sigma_c = reshape(sqrt(Pc(i,i,:)), 1, []);
            plot(tdata, -2*sigma_c, 'LineStyle', '--', ...
                'Color', 'green');
            plot(tdata, 2*sigma_c, 'LineStyle', '--', ...
                'Color', 'green');
            percent_bound(i,2) = sum(abs(errorsi) < 2*sigma_c)/length(tdata);
        end
        plot(tdata, 2*sigma, 'LineStyle', '--', ...
            'Color', 'red');
        
        ylabel(sprintf("%s Error [%s]", components(i), units(i)))
    
        if scale2err && max(3*sigma) > max(abs(errorsi))*3
            ylim([-mean(abs(errorsi))*10, mean(abs(errorsi))*10])
        end
    end
    xlabel("Time [s]")
    if ~isempty(Pc)
        lgd = legend(["Filter Error", "Filter 2\sigma Bounds", "Consider 2\sigma Bounds"], 'NumColumns', 3);
    else
        lgd = legend(["State Error", "Consider 2\sigma Bounds"], 'NumColumns', 2);
    end
    lgd.Layout.Tile = 'north';
    sgtitle(title_str);
end