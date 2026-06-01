function create_combined_window_FR_comparison(data, sess_types, neuron_types, ageCols, figDir)
% CREATE_COMBINED_WINDOW_FR_COMPARISON  Mean dFR (mono+poly averaged) per
% region (PL, IL), Y vs O, by neuron type. Sessions pooled. Per-rat means
% with jittered scatter.
% SS 2026

regions       = {'PL', 'IL'};
region_labels = {'Prelimbic', 'Infralimbic'};
ages          = {'Y', 'O'};

x_pos.INT = [1.0, 1.75];
x_pos.PYR = [3.25, 4.0];
bw = 0.65;

all_rm = [];
for ri = 1:2
    for ti = 1:numel(neuron_types)
        for ai = 1:2
            rm = collect_rat_means(data, sess_types, ages{ai}, neuron_types{ti}, regions{ri});
            all_rm = [all_rm; rm(:)];
        end
    end
end
ylim_u = get_unified_ylim({all_rm}, 0.25);

figure('Position', [100, 100, 800, 500]);
for ri = 1:2
    region = regions{ri};
    subplot(1, 2, ri); hold on;

    for ti = 1:numel(neuron_types)
        type = neuron_types{ti};
        xp = x_pos.(type);
        for ai = 1:2
            a = ages{ai};
            rm = collect_rat_means(data, sess_types, a, type, region);
            if isempty(rm); continue; end
            x = xp(ai);
            plot_scatter_with_jitter(gca, rm, x, ageCols.(a), ...
                'PRM_size', 80, 'PRM_jitter', 0.08, 'PRM_alpha', 0.7, 'PRM_edge', 'k');
            plot_bar_with_sem(gca, rm, x, bw, ageCols.(a), ...
                'PRM_face_alpha', 0.4, 'PRM_edge', 'k', 'PRM_eb_color', 'k');
        end
    end

    plot([2.625, 2.625], ylim_u, 'k--', 'LineWidth', 1);
    xlim([0.5, 4.75]); ylim(ylim_u);
    set(gca, 'XTick', [1.375, 3.625], 'XTickLabel', neuron_types);
    ylabel('\Delta FR (Hz)');
    title(region_labels{ri});
    box off;
    if exist('pubify_figure_axis_robust', 'file'); pubify_figure_axis_robust(16, 16); end

    if ri == 1
        h1 = plot(NaN, NaN, 's', 'MarkerFaceColor', ageCols.Y, 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
        h2 = plot(NaN, NaN, 's', 'MarkerFaceColor', ageCols.O, 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
        legend([h1, h2], {'Y', 'O'}, 'Location', 'northeast');
    end
end

sgtitle('Mean \Delta FR (Mono + Poly)');
saveas(gcf, fullfile(figDir, 'combined_window_FR_by_age_type.png'));
close(gcf);
end


function rm = collect_rat_means(data, sess_types, a, type, region)
rats_all = {}; vals_all = [];
for si = 1:numel(sess_types)
    d = data.(sess_types{si}).(a).(type);
    if isempty(d.rat); continue; end
    rmask = strcmp(d.region, region);
    if ~any(rmask); continue; end
    rats_all = [rats_all; d.rat(rmask)];
    vals_all = [vals_all; (d.win1(rmask) + d.win2(rmask))/2];
end
if isempty(vals_all); rm = []; return; end
rm = calculate_rat_means(rats_all, vals_all);
end
