function plot_iHC_vs_vHC_with_means(data, COLORS)
% PLOT_IHC_VS_VHC_WITH_MEANS  Within-rat iHC vs vHC comparison, individual
% lines transparent + mean +/- SEM bold. Combined analysis (both L2 and L5
% present) averages L2 and L5 per rat first. Creates two figures
% (Amplitude, Latency), each 1x4.
% SS 2026

if isfield(data, 'layers') && any(strcmp(data.layers, 'L2')) && any(strcmp(data.layers, 'L5'))
    data = average_layers_per_rat(data);
end

regions = {'PL', 'IL'};
windows = {'Mono', 'Poly'};

amp_all = []; lat_all = [];
for r = 1:2
    for w = 1:2
        amp_all = [amp_all; data.amplitude(:, w, r, 1); data.amplitude(:, w, r, 2)];
        lat_all = [lat_all; data.latency(:, w, r, 1);   data.latency(:, w, r, 2)];
    end
end
amp_ylim = pad_ylim(amp_all);
lat_ylim = pad_ylim(lat_all);

draw_figure(data, COLORS, regions, windows, 'amplitude', amp_ylim, 'Amplitude (mV)', ...
    'Within-Rat: iHC vs vHC - Amplitude (Negative Deflections)');
draw_figure(data, COLORS, regions, windows, 'latency', lat_ylim, 'Time to Peak (ms)', ...
    'Within-Rat: iHC vs vHC - Latency (Negative Deflections)');
end


function draw_figure(data, COLORS, regions, windows, field, ylim_vec, ylab, sg)
figure('Position', [100 100 1200 300], 'Name', ylab);
panel = 1;
for r = 1:2
    for w = 1:2
        subplot(1, 4, panel);
        d1 = data.(field)(:, w, r, 1);
        d2 = data.(field)(:, w, r, 2);
        plot_paired_comparison_panel(d1, d2, data.young_idx, data.old_idx, COLORS, ...
            ylim_vec, sprintf('%s - %ssynaptic', regions{r}, windows{w}), ...
            ylab, {'iHC', 'vHC'}, panel == 1);
        panel = panel + 1;
    end
end
sgtitle(sg);
end
