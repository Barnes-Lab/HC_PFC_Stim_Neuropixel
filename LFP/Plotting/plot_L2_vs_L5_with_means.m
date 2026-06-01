function plot_L2_vs_L5_with_means(data, COLORS)
% PLOT_L2_VS_L5_WITH_MEANS  Within-rat L2 vs L5 comparison (Combined
% analysis only: each rat must have both L2 and L5 rows). Two figures
% (Amplitude, Latency), each 1x4 (region x window). Sessions averaged.
% SS 2026

if ~isfield(data, 'layers')
    warning('Layer information not available - skipping L2 vs L5 plot'); return;
end
if ~any(strcmp(data.layers, 'L2')) || ~any(strcmp(data.layers, 'L5'))
    warning('Need both L2 and L5 data - skipping'); return;
end

regions = {'PL', 'IL'};
windows = {'Mono', 'Poly'};

unique_rats = unique(data.rat_ids);
n_rats      = numel(unique_rats);
rat_ages    = cell(n_rats, 1);
for i = 1:n_rats
    rat_ages{i} = data.ages{find(strcmp(data.rat_ids, unique_rats{i}), 1)};
end
young_idx = find(strcmp(rat_ages, 'Y'));
old_idx   = find(strcmp(rat_ages, 'O'));

% Collect for shared y-limits
amp_all = []; lat_all = [];
for r = 1:2
    for w = 1:2
        [L2a, L5a] = build_layer_arrays(data, unique_rats, n_rats, 'amplitude', w, r);
        [L2l, L5l] = build_layer_arrays(data, unique_rats, n_rats, 'latency',   w, r);
        amp_all = [amp_all; L2a; L5a];
        lat_all = [lat_all; L2l; L5l];
    end
end
amp_ylim = pad_ylim(amp_all);
lat_ylim = pad_ylim(lat_all);

draw_fig(data, unique_rats, n_rats, young_idx, old_idx, COLORS, regions, windows, ...
    'amplitude', amp_ylim, 'Amplitude (mV)', ...
    'Within-Rat: L2 vs L5 - Amplitude (Negative Deflections)');
draw_fig(data, unique_rats, n_rats, young_idx, old_idx, COLORS, regions, windows, ...
    'latency', lat_ylim, 'Time to Peak (ms)', ...
    'Within-Rat: L2 vs L5 - Latency (Negative Deflections)');
end


function draw_fig(data, unique_rats, n_rats, young_idx, old_idx, COLORS, regions, windows, ...
    field, ylim_vec, ylab, sg)
figure('Position', [100 100 1200 300], 'Name', ylab);
panel = 1;
for r = 1:2
    for w = 1:2
        subplot(1, 4, panel);
        [L2, L5] = build_layer_arrays(data, unique_rats, n_rats, field, w, r);
        plot_paired_comparison_panel(L2, L5, young_idx, old_idx, COLORS, ylim_vec, ...
            sprintf('%s - %ssynaptic', regions{r}, windows{w}), ylab, ...
            {'L2', 'L5'}, panel == 1);
        panel = panel + 1;
    end
end
sgtitle(sg);
end


function [L2_arr, L5_arr] = build_layer_arrays(data, unique_rats, n_rats, field, w, r)
% Per-rat L2 and L5 values for one (window, region), session-averaged.
L2_arr = nan(n_rats, 1);
L5_arr = nan(n_rats, 1);
for i = 1:n_rats
    L2_idx = find(strcmp(data.rat_ids, unique_rats{i}) & strcmp(data.layers, 'L2'), 1);
    L5_idx = find(strcmp(data.rat_ids, unique_rats{i}) & strcmp(data.layers, 'L5'), 1);
    if ~isempty(L2_idx) && ~isempty(L5_idx)
        L2_arr(i) = nanmean(data.(field)(L2_idx, w, r, :));
        L5_arr(i) = nanmean(data.(field)(L5_idx, w, r, :));
    end
end
end
