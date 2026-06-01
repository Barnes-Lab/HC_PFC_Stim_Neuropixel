function create_classification_summary_plot(classifications, all_SP_filtered, figDir, apply_artifact_rejection)
%Summary of neuron classification: 
% 1) scatter : trough-to-peak vs half-width (and pie chart with
% unclassifiedd)
% 2) mean waveforms by type 
% 3) firing-rate histograms 
% SS 2026

positive_peak_threshold = 0.5;
if nargin < 4; apply_artifact_rejection = false; end

classColors.INT          = [0.8, 0.2, 0.2];
classColors.PYR          = [0, 0.447, 0.741];
classColors.UNCLASSIFIED = [0.5, 0.5, 0.5];

if apply_artifact_rejection
    plot_idx = classifications.artifact_free;
else
    plot_idx = true(length(classifications.type), 1);
end

%% Scatter + pie chart
figure('Position', [100, 100, 600, 900]);
subplot(1, 2, 1); hold on;
for type = {'INT', 'PYR', 'UNCLASSIFIED'}
    idx = strcmp(classifications.type, type{1}) & plot_idx;
    if sum(idx) > 0
        scatter(classifications.half_width(idx), classifications.trough_to_peak(idx), ...
            40, classColors.(type{1}), 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', type{1});
    end
end
xlabel('Half-width (ms)'); ylabel('Trough-to-peak (ms)');
title('Waveform Shape Distribution');
legend('Location', 'best');
xlim([0, 0.8]); ylim([0, 1.2]); box off;
pubify_figure_axis_robust(20, 20);

subplot(1, 2, 2);
INT_count = sum(strcmp(classifications.type, 'INT')          & plot_idx);
PYR_count = sum(strcmp(classifications.type, 'PYR')          & plot_idx);
UNC_count = sum(strcmp(classifications.type, 'UNCLASSIFIED') & plot_idx);
counts    = [INT_count, PYR_count, UNC_count];
labels    = {'INT', 'PYR', 'UNCLASSIFIED'};
colors    = [classColors.INT; classColors.PYR; classColors.UNCLASSIFIED];

pie(counts);
colormap(gca, colors);
title_str = sprintf('Waveform Classification (n=%d)', sum(counts));
if apply_artifact_rejection; title_str = [title_str, ' - Artifact-free']; end
ht = title(title_str);
ht.Position = [ht.Position(1), ht.Position(2) + 0.1*ht.Position(2), ht.Position(3)];
legend(arrayfun(@(i) sprintf('%s: %d (%.1f%%)', labels{i}, counts(i), 100*counts(i)/sum(counts)), ...
    1:length(labels), 'UniformOutput', false), 'Location', 'southoutside', 'FontSize', 16);
pubify_figure_axis_robust(16, 16);

sgtitle('Neuron Classification Summary');
pubify_figure_axis_robust(16, 16);
saveas(gcf, fullfile(figDir, 'classification_summary.png'));
saveas(gcf, fullfile(figDir, 'classification_summary.fig'));
close(gcf);

plot_mean_waveforms_by_type(classifications, all_SP_filtered, figDir, classColors, ...
    apply_artifact_rejection, positive_peak_threshold);
plot_firing_rate_histograms(classifications, all_SP_filtered, figDir, classColors, ...
    apply_artifact_rejection);
end











