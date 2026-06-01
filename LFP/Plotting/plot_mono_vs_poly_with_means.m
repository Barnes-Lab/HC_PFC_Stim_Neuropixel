function plot_mono_vs_poly_with_means(data, COLORS)
%  Within-rat Mono vs Poly comparison 
% All Channels, PL, IL). 
% Combined analysis averages L2 and L5 per rat first.
% Sessions averaged
% SS 2026

if isfield(data, 'layers') && any(strcmp(data.layers, 'L2')) && any(strcmp(data.layers, 'L5'))
    data = average_layers_per_rat(data);
end

conditions = {'All Channels', 'PL', 'IL'};

all_data = [];
for cond = 1:3
    [m, p] = get_mono_poly(data, cond);
    all_data = [all_data; m; p];
end
y_lim = pad_ylim(all_data, 0.10, 0.10);

for cond = 1:3
    subplot(1, 3, cond);
    [mono_data, poly_data] = get_mono_poly(data, cond);
    plot_paired_comparison_panel(mono_data, poly_data, data.young_idx, data.old_idx, ...
        COLORS, y_lim, conditions{cond}, 'Amplitude (mV)', {'Mono', 'Poly'}, cond == 1);
end
sgtitle('Within-Rat: Monosynaptic vs Polysynaptic (Negative Deflections)');
end


function [mono_d, poly_d] = get_mono_poly(data, cond)
% Region/session-averaged mono and poly amplitude per rat for one condition.
switch cond
    case 1
        mono_d = squeeze(nanmean(data.amplitude(:, 1, :, :), [3 4]));
        poly_d = squeeze(nanmean(data.amplitude(:, 2, :, :), [3 4]));
    case 2
        mono_d = squeeze(nanmean(data.amplitude(:, 1, 1, :), 4));
        poly_d = squeeze(nanmean(data.amplitude(:, 2, 1, :), 4));
    case 3
        mono_d = squeeze(nanmean(data.amplitude(:, 1, 2, :), 4));
        poly_d = squeeze(nanmean(data.amplitude(:, 2, 2, :), 4));
end
end
