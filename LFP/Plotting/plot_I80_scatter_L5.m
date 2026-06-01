function plot_I80_scatter_L5(I80_TABLE, PARAMS, output_dir, included_rats)
% Scatter of I80 (uA) for Young vs Old, split into iHC and vHC subplots.
% Layout matches the original plot_I80_scatter_L5 from RUN_LF_IO_CURVES;
% only change is that included_rats is now a parameter rather than
% hardcoded inside. Reads I80 from the chosen-shank table built by
% the new I80 pipeline (mean_amplitude fits).
%
% Two-sample t-test per stim site is returned
%
%
% SS 2026

if nargin < 4 || isempty(included_rats)
    included_rats = unique(I80_TABLE.Rat);
end

rat_mask = ismember(I80_TABLE.Rat, included_rats);
T = I80_TABLE(rat_mask, :);

fig = figure('Position', [100, 100, 1000, 500]);
stim_sites = {'iHC', 'vHC'};

for s = 1:2
    subplot(1, 2, s); hold on;
    site = stim_sites{s};

    rows = T(strcmp(T.Stim_Site, site), :);
    yI80 = rows.I80(strcmp(rows.Age, 'Y'));
    oI80 = rows.I80(strcmp(rows.Age, 'O'));

    jit = 0.1;
    yx = ones(size(yI80))   + (rand(size(yI80))   - 0.5) * jit;
    ox = 2 * ones(size(oI80)) + (rand(size(oI80)) - 0.5) * jit;

    scatter(yx, yI80, 180, PARAMS.colors.young, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);
    scatter(ox, oI80, 180, PARAMS.colors.old, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);

    if ~isempty(yI80)
        plot([0.7, 1.3], [mean(yI80), mean(yI80)], '-', ...
            'Color', PARAMS.colors.young, 'LineWidth', 4);
    end
    if ~isempty(oI80)
        plot([1.7, 2.3], [mean(oI80), mean(oI80)], '-', ...
            'Color', PARAMS.colors.old, 'LineWidth', 4);
    end

    xlim([0.5, 2.5]); ylim([100, 600]);
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'Young', 'Old'});
    ylabel('I80 (\muA)');
    title(site, 'FontWeight', 'bold');
    pubify_figure_axis_robust(20, 20);
    box off;

    if numel(yI80) >= 2 && numel(oI80) >= 2
        [~, p, ~, stats] = ttest2(yI80, oI80);
        fprintf('%s: Y N=%d M=%.1f SD=%.1f | O N=%d M=%.1f SD=%.1f | t(%d)=%.3f p=%.4f\n', ...
            site, numel(yI80), mean(yI80), std(yI80), ...
            numel(oI80), mean(oI80), std(oI80), stats.df, stats.tstat, p);
    end
end

sgtitle('Stim Intensity at 80% of Monosynaptic Plateau Amplitude', ...
    'FontSize', 20, 'FontWeight', 'bold');

base = 'I80_Scatter_L5';
saveas(fig, fullfile(output_dir, [base '.png']));
close(fig);
fprintf('Saved %s.png\n', base);

end
