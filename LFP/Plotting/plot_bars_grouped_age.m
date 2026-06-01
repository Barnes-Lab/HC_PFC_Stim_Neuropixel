function plot_bars_grouped_age(data, COLORS, save_dir)
% Four grouped bar plots with Young/Old for IL and PL
% at each region (PL, IL):
% ono Amplitude, Mono Time to Peak, Poly Amplitude, Poly Time to Peak.
% 
% SS 2026

if isfield(data, 'layers') && any(strcmp(data.layers, 'L2')) && any(strcmp(data.layers, 'L5'))
    data = average_layers_per_rat(data);
end

regions = {'PL', 'IL'};

amp_all = []; lat_all = [];
for window = 1:2
    for r = 1:2
        amp_all = [amp_all;
            squeeze(nanmean(data.amplitude(data.young_idx, window, r, :), 4));
            squeeze(nanmean(data.amplitude(data.old_idx,   window, r, :), 4))];
        lat_all = [lat_all;
            squeeze(nanmean(data.latency(data.young_idx, window, r, :), 4));
            squeeze(nanmean(data.latency(data.old_idx,   window, r, :), 4))];
    end
end
amp_ylim = pad_ylim(amp_all, 0.25, 0.35);
lat_ylim = pad_ylim(lat_all, 0.25, 0.35);

win_names    = {'MONO', 'POLY'};
win_titles   = {'Monosynaptic', 'Polysynaptic'};
metric_data  = {data.amplitude, data.latency};
metric_ylim  = {amp_ylim, lat_ylim};
metric_label = {'Amplitude (mV)', 'Time to Peak (ms)'};
metric_name  = {'Amplitude', 'TimeToPeak'};

bar_width   = 0.28;
x_positions = [1, 2];

for window = 1:2
    for metric = 1:2
        figure('Position', [100 100 300 400]); hold on;
        for r = 1:2
            young_d = squeeze(nanmean(metric_data{metric}(data.young_idx, window, r, :), 4));
            old_d   = squeeze(nanmean(metric_data{metric}(data.old_idx,   window, r, :), 4));
            draw_one_region(x_positions(r), bar_width, young_d, old_d, COLORS);
        end
        xlim([0.5, 2.5]); ylim(metric_ylim{metric});
        set(gca, 'XTick', x_positions, 'XTickLabel', regions);
        ylabel(metric_label{metric});
        title(sprintf('%s %s (Averaged Across Sessions)', win_titles{window}, metric_name{metric}));
        legend({'Young', 'Old'}, 'Location', 'eastoutside');
        pubify_figure_axis_robust(20, 20);

        fname = sprintf('Figure_%s_%s', win_names{window}, metric_name{metric});
        saveas(gcf, fullfile(save_dir, [fname '.png']));
        saveas(gcf, fullfile(save_dir, [fname '.fig']));
        close(gcf);
    end
end
fprintf('Grouped bar plots complete (4 figures)\n');
end


function draw_one_region(x_center, bar_width, young_d, old_d, COLORS)
% One PL or IL region: scatter both ages first, then bars + errorbars on top.
ax = gca;
xY = x_center - bar_width/2;
xO = x_center + bar_width/2;

plot_scatter_with_jitter(ax, young_d(~isnan(young_d)), xY, COLORS.young, ...
    'PRM_MARKER_SIZE', 80, 'PRM_JITTER', 0.08, 'PRM_LIGHTEN', false, ...
    'PRM_EDGE', 'k', 'PRM_FACE_ALPHA', 0.7);
plot_scatter_with_jitter(ax, old_d(~isnan(old_d)),   xO, COLORS.old, ...
    'PRM_MARKER_SIZE', 80, 'PRM_JITTER', 0.08, 'PRM_LIGHTEN', false, ...
    'PRM_EDGE', 'k', 'PRM_FACE_ALPHA', 0.7);

plot_bar_with_sem(ax, young_d, xY, COLORS.young, ...
    'PRM_WIDTH', bar_width, 'PRM_FACE_ALPHA', 0.4, 'PRM_EDGE', 'k', ...
    'PRM_LINE_WIDTH', 1.5, 'PRM_ERR_LW', 1.5);
plot_bar_with_sem(ax, old_d, xO, COLORS.old, ...
    'PRM_WIDTH', bar_width, 'PRM_FACE_ALPHA', 0.4, 'PRM_EDGE', 'k', ...
    'PRM_LINE_WIDTH', 1.5, 'PRM_ERR_LW', 1.5);
end
