function plot_psth_depth_heatmap(timeBins_ms, depthBins, allP, titleStr, PL_depth)
%Imagesc of depth x time PSTH with stim line at
% t=0, PL/IL boundary at PL_depth, white-to-red colormap, colorbar in Hz.


imagesc(timeBins_ms, depthBins(1:end-1), allP);
set(gca, 'YDir', 'reverse');
hold on;

nD = length(depthBins) - 1;
plot(zeros(1, nD), depthBins(1:end-1), 'Color', '#A2142F', 'LineWidth', 3);
colormap(gca, cmap_white_to_red(256, 1));

xlabel('Time from stim (ms)');
xlim([min(timeBins_ms), max(timeBins_ms)]);
ylabel('Depth (\mum)');
title(titleStr);
box off;

h = colorbar;
h.Label.String = '\Delta FR (Hz)';

plot([min(timeBins_ms), max(timeBins_ms)], [PL_depth, PL_depth], 'k:', 'LineWidth', 3);
set(gca, 'FontSize', 14);
end
