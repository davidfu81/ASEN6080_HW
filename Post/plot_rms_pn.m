function figure_rms_pn = plot_rms_pn(sigmas_m, rms_vals, ylabels, title_str)
    figure_rms_pn = figure;
    tiledlayout(size(rms_vals, 1),1)

    for i = 1:size(rms_vals,1)
        nexttile; grid on;
        plot(sigmas_m, rms_vals(i,:))
        xscale('log')
        ylabel(ylabels(i))
    end
    
    xlabel("Process Noise Sigma [m/s^2]")
    
    sgtitle(title_str)
end