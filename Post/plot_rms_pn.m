function figure_rms_pn = plot_rms_pn(sigmas_m, rms_vals, ylabels, title_str)
    figure_rms_pn = figure;
    set(gcf, "Position", [100 100 750 500])
    tl = tiledlayout(size(rms_vals, 1),1);
    tl.TileSpacing = 'compact';
    tl.Padding = 'loose';
    for i = 1:size(rms_vals,1)
        nexttile; grid on;
        plot(sigmas_m, rms_vals(i,:))
        xscale('log')
        ylabel(ylabels(i))
        % if max(rms_vals(i,:)) > 10*mean(rms_vals(i,:))
        %     ylim([0, mean(rms_vals(i,:))*2])
        % end
        yscale('log')
        xlim([min(sigmas_m), max(sigmas_m)])
    end
    
    xlabel("Process Noise Sigma [m/s^2]")
    
    sgtitle(title_str)
end