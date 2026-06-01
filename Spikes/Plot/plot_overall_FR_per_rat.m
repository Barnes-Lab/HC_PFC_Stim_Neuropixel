function plot_overall_FR_per_rat(data, sess_types, mode, ageCols, figDir, recruit_thresh)
% PLOT_OVERALL_FR_PER_RAT  Per-rat combined-window dFR (PYR), Y
%   Rows:PL, IL, PL+IL combined
%   Columns: Combined sess (iHC+vHC), iHC only, vHC only
%
% 
% Modes:
%   'all'  - all neurons contribute to each rat's mean
%   'responders'- only neurons > thresh per window included
%
% Combined window per neuron = [mono_start, poly_end] 
%
% SS 2026

regions       = {'PL', 'IL', 'all'};
region_titles = {'PL', 'IL', 'PL+IL'};
sess_filters  = {'all', 'iHC', 'vHC'};
sess_titles   = {'Combined', 'iHC', 'vHC'};

resp_only = strcmp(mode, 'responders');
rats_Y = unique_rats_for_age(data, sess_types, 'PYR', 'Y');
rats_O = unique_rats_for_age(data, sess_types, 'PYR', 'O');

figure('Position', [100, 100, 1100, 900]);
for ri = 1:3
    for sj = 1:3
        ax = subplot(3, 3, (ri - 1) * 3 + sj);
        vY = aggregate_win_combined_per_rat(data, sess_types, 'Y', ...
            sess_filters{sj}, regions{ri}, resp_only, recruit_thresh, rats_Y);
        vO = aggregate_win_combined_per_rat(data, sess_types, 'O', ...
            sess_filters{sj}, regions{ri}, resp_only, recruit_thresh, rats_O);
        vY = vY(~isnan(vY)); vO = vO(~isnan(vO));

        draw_age_panel(ax, vY, vO, ageCols);
        ylabel(ax, '\Delta FR (Hz)');
        title(ax, sprintf('%s - %s (n_Y=%d, n_O=%d)', ...
            region_titles{ri}, sess_titles{sj}, numel(vY), numel(vO)));
    end
end

mode_str = mode; if resp_only; mode_str = 'responders only'; end
sgtitle(sprintf('Overall combined-window \\Delta FR per rat [%s, PYR]', mode_str));

base = fullfile(figDir, sprintf('overall_FR_per_rat_%s_PYR', mode));
saveas(gcf, [base '.png']);
saveas(gcf, [base '.fig']);
close(gcf);
end


function draw_age_panel(ax, vY, vO, ageCols)
axes(ax); hold(ax, 'on');
x_Y = 1; x_O = 2;

plot_scatter_with_jitter(ax, vY, x_Y, ageCols.Y, 'PRM_MARKER_SIZE', 55, 'PRM_JITTER', 0.18);
plot_scatter_with_jitter(ax, vO, x_O, ageCols.O, 'PRM_MARKER_SIZE', 55, 'PRM_JITTER', 0.18);

plot_bar_with_sem(ax, vY, x_Y, ageCols.Y);
plot_bar_with_sem(ax, vO, x_O, ageCols.O);

xlim(ax, [0.4, 2.6]);
set(ax, 'XTick', [x_Y, x_O], 'XTickLabel', {'Y', 'O'});
box(ax, 'off');
if exist('pubify_figure_axis_robust', 'file'); pubify_figure_axis_robust(12, 12); end
end
