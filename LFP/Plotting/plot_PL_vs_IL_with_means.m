function plot_PL_vs_IL_with_means(data, COLORS)
%  Within-rat PL vs IL comparison. 
% Combined analysis (both L2 and L5) averages L2 and L5 per rat first.
% (Amplitude, Latency)
% SS 2026

if isfield(data, 'layers') && any(strcmp(data.layers, 'L2')) && any(strcmp(data.layers, 'L5'))
    data = average_layers_per_rat(data);
end

sessions = {'iHC-Mono', 'iHC-Poly', 'vHC-Mono', 'vHC-Poly', 'All-Mono', 'All-Poly'};

[amp_all, lat_all] = collect_pl_il_values(data);
amp_ylim = pad_ylim(amp_all);
lat_ylim = pad_ylim(lat_all);

draw_pl_il_fig(data, COLORS, sessions, 'amplitude', amp_ylim, 'Amplitude (mV)', ...
    'Within-Rat: PL vs IL - Amplitude (Negative Deflections)');
draw_pl_il_fig(data, COLORS, sessions, 'latency', lat_ylim, 'Time to Peak (ms)', ...
    'Within-Rat: PL vs IL - Latency (Negative Deflections)');
end


function draw_pl_il_fig(data, COLORS, sessions, field, ylim_vec, ylab, sg)
figure('Position', [100 100 1800 300], 'Name', ylab);
panel = 1;
for sess = 1:3
    for w = 1:2
        subplot(1, 6, panel);
        if sess == 3
            PL_data = squeeze(nanmean(data.(field)(:, w, 1, :), 4));
            IL_data = squeeze(nanmean(data.(field)(:, w, 2, :), 4));
        else
            PL_data = data.(field)(:, w, 1, sess);
            IL_data = data.(field)(:, w, 2, sess);
        end
        plot_paired_comparison_panel(PL_data, IL_data, data.young_idx, data.old_idx, ...
            COLORS, ylim_vec, sessions{panel}, ylab, {'PL', 'IL'}, false);
        panel = panel + 1;
    end
end
sgtitle(sg);
end


function [amp_all, lat_all] = collect_pl_il_values(data)
amp_all = []; lat_all = [];
for sess = 1:3
    for w = 1:2
        if sess == 3
            amp_all = [amp_all; squeeze(nanmean(data.amplitude(:, w, 1, :), 4));
                                squeeze(nanmean(data.amplitude(:, w, 2, :), 4))];
            lat_all = [lat_all; squeeze(nanmean(data.latency(:, w, 1, :), 4));
                                squeeze(nanmean(data.latency(:, w, 2, :), 4))];
        else
            amp_all = [amp_all; data.amplitude(:, w, 1, sess); data.amplitude(:, w, 2, sess)];
            lat_all = [lat_all; data.latency(:, w, 1, sess);   data.latency(:, w, 2, sess)];
        end
    end
end
end
