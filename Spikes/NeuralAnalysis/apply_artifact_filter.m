function neuron_classifications = apply_artifact_filter(neuron_classifications, all_SP_filtered, figDir) %#ok<INUSD>
% APPLY_ARTIFACT_FILTER  Mark classified neurons whose mean waveform has a
% large positive deflection in a window around the average
% trough )+/-5ms. This is to check alignment to peak directly 
%
% INPUT
%   neuron_classifications - output of classify_neurons_filtered
%   all_SP_filtered      
%   figDir 
%
% OUTPUT
%   neuron_classification
%
% SS 2026

positive_peak_threshold = 0.95;
artifact_window         = 5;

classified_idx = strcmp(neuron_classifications.type, 'INT') | ...
                 strcmp(neuron_classifications.type, 'PYR');

trough_times = [];
for u = 1:length(all_SP_filtered)
    if classified_idx(u) && isfield(all_SP_filtered(u), 'WV') && ...
            isfield(all_SP_filtered(u).WV, 'mWV')
        wf = all_SP_filtered(u).WV.mWV(:, 1)';
        wf = wf - mean(wf(15:25));
        [~, trough_idx] = min(wf);
        trough_times(end+1) = trough_idx;
    end
end

if isempty(trough_times)
    fprintf('  Warning: No valid waveforms found\n'); return;
end
avg_trough_time = round(mean(trough_times));

artifact_count = 0;
for u = 1:length(all_SP_filtered)
    if classified_idx(u) && isfield(all_SP_filtered(u), 'WV') && ...
            isfield(all_SP_filtered(u).WV, 'mWV')
        wf = all_SP_filtered(u).WV.mWV(:, 1)';
        wf = wf - mean(wf(15:25));
        wf_length = length(wf);
        start_idx = max(1, avg_trough_time - artifact_window);
        end_idx   = min(wf_length, avg_trough_time + artifact_window);

        trough_val   = min(wf);
        has_artifact = any(wf(start_idx:end_idx) > positive_peak_threshold * abs(trough_val));
        if has_artifact
            neuron_classifications.artifact_free(u) = false;
            artifact_count = artifact_count + 1;
        end
    end
end

fprintf('  %d/%d classified neurons rejected as artifacts\n', ...
    artifact_count, sum(classified_idx));
end
