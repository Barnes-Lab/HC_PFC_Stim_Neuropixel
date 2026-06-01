function plot_mean_waveforms_by_type(classifications, all_SP_filtered, figDir, classColors, ...
    apply_artifact_rejection, positive_peak_threshold)
types_to_plot = {'INT', 'PYR'};
for ti = 1:length(types_to_plot)
    type = types_to_plot{ti};
    type_idx_logical = strcmp(classifications.type, type);
    if apply_artifact_rejection; type_idx_logical = type_idx_logical & classifications.artifact_free; end

    figure('Position', [100, 100, 800, 600]); hold on;
    [mean_wf, individual_wfs, ~] = extract_and_normalize_waveforms( ...
        all_SP_filtered, type_idx_logical, positive_peak_threshold);
    if ~isempty(individual_wfs)
        plot_waveform_overlay(individual_wfs, mean_wf, classColors.(type), type);
    end
    axis off; set(gca, 'XColor', 'none', 'YColor', 'none');
    saveas(gcf, fullfile(figDir, sprintf('mean_waveforms_%s.svg', type)));
    close(gcf);

    if apply_artifact_rejection
        rejected_idx = strcmp(classifications.type, type) & ~classifications.artifact_free;
        if sum(rejected_idx) > 0
            figure('Position', [100, 100, 800, 600]); hold on;
            [mean_wf_rej, individual_wfs_rej, ~] = extract_and_normalize_waveforms( ...
                all_SP_filtered, rejected_idx, Inf);
            if ~isempty(individual_wfs_rej)
                plot_waveform_overlay(individual_wfs_rej, mean_wf_rej, classColors.(type) * 0.6, type);
            end
            axis off; set(gca, 'XColor', 'none', 'YColor', 'none');
            saveas(gcf, fullfile(figDir, sprintf('mean_waveforms_%s_REJECTED.svg', type)));
            close(gcf);
            fprintf('  Saved %d rejected %s waveforms\n', sum(rejected_idx), type);
        end
    end
end
end
