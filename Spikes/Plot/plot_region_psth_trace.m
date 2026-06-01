function plot_region_psth_trace(timeBins_ms, region_psth, titleStr, color)
% Mean +/- SEM PSTH across depth bins for one


mean_psth = mean(region_psth, 1);
sem_psth  = std(region_psth, 0, 1) / sqrt(size(region_psth, 1));

hold on;
fill([timeBins_ms, fliplr(timeBins_ms)], ...
    [mean_psth + sem_psth, fliplr(mean_psth - sem_psth)], ...
    color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(timeBins_ms, mean_psth, 'Color', color, 'LineWidth', 2);
plot([0 0], ylim, 'r-', 'LineWidth', 2);
plot(xlim, [0 0], 'k:', 'LineWidth', 1);

xlabel('Time from stim (ms)');
ylabel('\Delta FR (Hz)');
title(titleStr);
xlim([timeBins_ms(1), timeBins_ms(end)]);
box off;
set(gca, 'FontSize', 12);
end
