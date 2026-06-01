function plot_firing_rate_histograms(classifications, all_SP_filtered, figDir, classColors, ...
    apply_artifact_rejection)
if apply_artifact_rejection; plot_idx = classifications.artifact_free;
else;                        plot_idx = true(length(classifications.type), 1);
end

fr_values = zeros(length(all_SP_filtered), 1);
has_fr    = false(length(all_SP_filtered), 1);
for u = 1:length(all_SP_filtered)
    if isfield(all_SP_filtered(u), 'fr') && ~isempty(all_SP_filtered(u).fr)
        fr_values(u) = all_SP_filtered(u).fr;
        has_fr(u)    = true;
    end
end

figure('Position', [100, 100, 800, 600]); hold on;
INT_idx = strcmp(classifications.type, 'INT') & has_fr & plot_idx;
PYR_idx = strcmp(classifications.type, 'PYR') & has_fr & plot_idx;
if sum(INT_idx) > 0 || sum(PYR_idx) > 0
    max_fr = max(fr_values(INT_idx | PYR_idx));
    bins   = 0:0.5:ceil(max_fr);
    if sum(INT_idx) > 0
        histogram(fr_values(INT_idx), bins, 'FaceColor', classColors.INT, ...
            'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'INT', 'Normalization', 'probability');
    end
    if sum(PYR_idx) > 0
        histogram(fr_values(PYR_idx), bins, 'FaceColor', classColors.PYR, ...
            'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'PYR', 'Normalization', 'probability');
    end
    xlabel('Firing Rate (Hz)'); ylabel('Proportion');
    legend('Location', 'best'); box off;
    pubify_figure_axis_robust(16, 16);
end
saveas(gcf, fullfile(figDir, 'firing_rate_histogram.png'));
close(gcf);
end
