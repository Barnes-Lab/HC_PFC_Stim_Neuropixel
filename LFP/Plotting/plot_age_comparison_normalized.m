function plot_age_comparison_normalized(I80_TABLE, MONO_DATA, PARAMS, ...
    output_dir, included_rats, amp_field)
%   Normalized IO curves by age group.
% Each rat is binned into [100,200,300,400,500,600] uA bins, then
% normalized by subtracting the rat's least-negative bin (baseline).
% Group means + SEM overlaid. iHC and vHC in side-by-side subplots.
%
% INPUT
%   I80_TABLE 
%   MONO_DATA- mono-detection output
%   PARAMS     -with .colors.young, .colors.old
%   output_dir  
%   included_rats
%   amp_field  (mean)
%
% SS 2026

if nargin < 5 || isempty(included_rats); included_rats = unique(I80_TABLE.Rat); end
if nargin < 6 || isempty(amp_field);     amp_field     = 'mean_amplitude';      end

T = I80_TABLE(ismember(I80_TABLE.Rat, included_rats), :);

stim_bins   = [100, 200, 300, 400, 500, 600];
bin_centers = [150, 250, 350, 450, 550];

fig = figure('Position', [100, 100, 1400, 600]);
stim_sites = {'iHC', 'vHC'};

for s = 1:2
    subplot(1, 2, s); hold on;
    site = stim_sites{s};
    rows = T(strcmp(T.Stim_Site, site), :);

    [yC, yA, ySE] = process_age_group(rows(strcmp(rows.Age, 'Y'), :), MONO_DATA, ...
        stim_bins, bin_centers, amp_field);
    [oC, oA, oSE] = process_age_group(rows(strcmp(rows.Age, 'O'), :), MONO_DATA, ...
        stim_bins, bin_centers, amp_field);

    draw_rat_traces(yC, bin_centers, PARAMS.colors.young);
    draw_rat_traces(oC, bin_centers, PARAMS.colors.old);
    draw_group_mean(yA, ySE, bin_centers, PARAMS.colors.young, 'o');
    draw_group_mean(oA, oSE, bin_centers, PARAMS.colors.old,   's');

    xlabel('Current (\muA)', 'FontSize', 12);
    ylabel('Normalized Amplitude (mV)', 'FontSize', 12);
    title(site, 'FontWeight', 'bold', 'FontSize', 13);
    xlim([80, 600]); ylim([-1.1, 0.1]);
    pubify_figure_axis_robust(20, 20);
    box off;

    if s == 1
        h = legend({'Young', 'Old'}, 'Location', 'southwest', 'Box', 'off');
        set(h, 'FontSize', 20);
    end
end

sgtitle('IL: Stimulus Response Normalized to Baseline (chosen shank)', ...
    'FontSize', 20, 'FontWeight', 'bold');
saveas(fig, fullfile(output_dir, 'Age_Comparison_Normalized_L5.png'));
close(fig);
fprintf('Saved Age_Comparison_Normalized_L5.png\n');
end


function draw_rat_traces(curves, bin_centers, col)
for i = 1:size(curves, 1)
    v = ~isnan(curves(i, :));
    if any(v)
        plot(bin_centers(v), curves(i, v), '-', 'Color', [col, 0.5], ...
            'LineWidth', 2, 'HandleVisibility', 'off');
    end
end
end


function draw_group_mean(avg, sem, bin_centers, col, marker)
v = ~isnan(avg);
if ~any(v); return; end
plot(bin_centers(v), avg(v), '-', 'Color', col, 'LineWidth', 5, 'HandleVisibility', 'off');
errorbar(bin_centers(v), avg(v), sem(v), marker, 'Color', col, 'LineWidth', 2, ...
    'MarkerSize', 12, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
    'CapSize', 10, 'LineStyle', 'none');
end
