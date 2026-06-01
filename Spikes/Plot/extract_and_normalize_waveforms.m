function [mean_wf, individual_wfs, valid_idx] = extract_and_normalize_waveforms( ...
    all_SP_filtered, idx, threshold)
% Extract waveforms, reject artifacts, normalize, upsample, align to mean trough.
individual_wfs = [];
valid_idx      = false(length(all_SP_filtered), 1);
fs             = 30000;
n_up_pts       = 1000;

trough_times = [];
for u = 1:length(all_SP_filtered)
    if idx(u) && isfield(all_SP_filtered(u), 'WV') && isfield(all_SP_filtered(u).WV, 'mWV')
        wf = all_SP_filtered(u).WV.mWV(:, 1)';
        wf = wf - mean(wf(15:25));
        [~, trough_idx] = min(wf);
        trough_times(end+1) = trough_idx;
    end
end
if isempty(trough_times); mean_wf = []; return; end
avg_trough_time = round(mean(trough_times));

artifact_window = 5;
first_wf_idx    = find(idx, 1);
if isempty(first_wf_idx); mean_wf = []; return; end
wf_length = length(all_SP_filtered(first_wf_idx).WV.mWV(:, 1)');
start_idx = max(1, avg_trough_time - artifact_window);
end_idx   = min(wf_length, avg_trough_time + artifact_window);

temp_waveforms      = [];
temp_neuron_indices = [];
for u = 1:length(all_SP_filtered)
    if idx(u) && isfield(all_SP_filtered(u), 'WV') && isfield(all_SP_filtered(u).WV, 'mWV')
        wf = all_SP_filtered(u).WV.mWV(:, 1)';
        wf = wf - mean(wf(15:25));
        trough_val   = min(wf);
        has_artifact = any(wf(start_idx:end_idx) > threshold * abs(trough_val));
        if ~has_artifact
            temp_waveforms      = [temp_waveforms; wf];
            temp_neuron_indices = [temp_neuron_indices; u];
        end
    end
end
if isempty(temp_waveforms)
    mean_wf = []; fprintf('  Warning: No waveforms passed artifact rejection\n'); return;
end

individual_wfs  = zeros(size(temp_waveforms, 1), n_up_pts);
trough_indices  = zeros(size(temp_waveforms, 1), 1);
tmp             = (1:wf_length) - 1;
wv_x_ms         = 1000 * (tmp / fs);
wv_x_ms_up      = linspace(wv_x_ms(1), wv_x_ms(end), n_up_pts);

for i = 1:size(temp_waveforms, 1)
    wave_up = interp1(wv_x_ms, temp_waveforms(i, :), wv_x_ms_up, 'spline');
    [~, trough_idx]   = min(wave_up);
    trough_indices(i) = trough_idx;
    individual_wfs(i, :) = wave_up / abs(min(wave_up));
end

mean_trough_idx = round(mean(trough_indices));
aligned_wfs     = zeros(size(individual_wfs));
for i = 1:size(individual_wfs, 1)
    shift = mean_trough_idx - trough_indices(i);
    if shift > 0
        aligned_wfs(i, (shift+1):end) = individual_wfs(i, 1:(end-shift));
        aligned_wfs(i, 1:shift)       = individual_wfs(i, 1);
    elseif shift < 0
        aligned_wfs(i, 1:(end+shift))   = individual_wfs(i, (1-shift):end);
        aligned_wfs(i, (end+shift+1):end) = individual_wfs(i, end);
    else
        aligned_wfs(i, :) = individual_wfs(i, :);
    end
end
individual_wfs = aligned_wfs;

for i = 1:length(temp_neuron_indices)
    valid_idx(temp_neuron_indices(i)) = true;
end

mean_wf = mean(individual_wfs, 1);
fprintf('  %d/%d waveforms included\n', size(individual_wfs, 1), sum(idx));
end
