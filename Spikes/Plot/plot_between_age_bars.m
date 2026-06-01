function plot_between_age_bars(data, sess_types, granularity, metric, ageCols, figDir, y_label, recruit_thresh)
% PLOT_BETWEEN_AGE_BARS  Between-age (Y vs O) bar plots, PYR-only.
% 7 conditions per panel:
%   Per-region group (4): PL-mono, PL-poly, IL-mono, IL-poly
%   Overall group   (3): Mono-all, Poly-all, All-all  (regions pooled)
% Vertical dashed line separates the two groups. Matches rows of the
% between-age post-hoc CSV.
%
% Layering: rat scatter (both ages) first, then bars + errorbars on top.
%
% Granularity:
%   'coarse'      subplot(1,1) sessions collapsed
%   'by_session'  subplot(1,2) iHC | vHC
%
% SS 2026

[panel_titles, sess_filters, layout, fig_size] = get_layout(granularity);

cond_specs = {
    'PL-mono',  'PL',  'mono', 1;
    'PL-poly',  'PL',  'poly', 1;
    'IL-mono',  'IL',  'mono', 1;
    'IL-poly',  'IL',  'poly', 1;
    'Mono-all', 'avg', 'mono', 2;
    'Poly-all', 'avg', 'poly', 2;
    'All-all',  'avg', 'avg',  2;
};
cond_labels = cond_specs(:, 1)';

rats_Y = unique_rats_for_age(data, sess_types, 'PYR', 'Y');
rats_O = unique_rats_for_age(data, sess_types, 'PYR', 'O');

figure('Position', [100, 100, fig_size(1), fig_size(2)]);
for p = 1:numel(sess_filters)
    sess = sess_filters{p};
    nC = size(cond_specs, 1);
    Y_vals = cell(1, nC); O_vals = cell(1, nC);
    for c = 1:nC
        region = cond_specs{c, 2};
        win    = cond_specs{c, 3};
        vY = get_per_rat_value(data, sess_types, 'Y', sess, region, win, metric, recruit_thresh, rats_Y);
        vO = get_per_rat_value(data, sess_types, 'O', sess, region, win, metric, recruit_thresh, rats_O);
        Y_vals{c} = vY(~isnan(vY));
        O_vals{c} = vO(~isnan(vO));
    end
    ax = subplot(layout(1), layout(2), p);
    draw_age_bars_panel(ax, Y_vals, O_vals, cond_labels, cell2mat(cond_specs(:, 4)'), ageCols);
    ylabel(ax, y_label);
    title(ax, panel_titles{p});
end

sgtitle(sprintf('Y vs O %s: %s [PYR]', strrep(granularity, '_', ' '), metric), 'Interpreter', 'none');
base = fullfile(figDir, sprintf('between_age_bars_%s_%s_PYR', granularity, metric));
saveas(gcf, [base '.png']);
saveas(gcf, [base '.fig']);
close(gcf);
end


function draw_age_bars_panel(ax, Y_vals, O_vals, cond_labels, group_ids, ageCols)
axes(ax); hold(ax, 'on');
nC    = numel(cond_labels);
x     = 1:nC;
width = 0.35;

for c = 1:nC; plot_scatter_with_jitter(ax, Y_vals{c}, x(c) - width/2, ageCols.Y); end
for c = 1:nC; plot_scatter_with_jitter(ax, O_vals{c}, x(c) + width/2, ageCols.O); end

for c = 1:nC; plot_bar_with_sem(ax, Y_vals{c}, x(c) - width/2, ageCols.Y, 'PRM_WIDTH', width); end
for c = 1:nC; plot_bar_with_sem(ax, O_vals{c}, x(c) + width/2, ageCols.O, 'PRM_WIDTH', width); end

for b = find(diff(group_ids) ~= 0)
    xline(ax, b + 0.5, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

xlim(ax, [0.5, nC + 0.5]);
set(ax, 'XTick', x, 'XTickLabel', cond_labels);
ax.XTickLabelRotation = 30;
box(ax, 'off');
if exist('pubify_figure_axis_robust', 'file'); pubify_figure_axis_robust(13, 13); end
end


function [titles, sess_filters, layout, fig_size] = get_layout(granularity)
switch granularity
    case 'coarse'
        titles = {'Sessions collapsed'}; sess_filters = {'all'};
        layout = [1, 1]; fig_size = [850, 550];
    case 'by_session'
        titles = {'iHC', 'vHC'};          sess_filters = {'iHC', 'vHC'};
        layout = [1, 2]; fig_size = [1500, 550];
    otherwise
        error('Unknown granularity: %s (use ''coarse'' or ''by_session'')', granularity);
end
end
