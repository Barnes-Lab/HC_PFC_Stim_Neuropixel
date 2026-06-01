function reg = run_region_detection_replace(Vb_all, T_VEC, region_chans, ...
    stim_idx, extras_pool, allDepths, t_b_start, t_b_end, ...
    mono_win, poly_win, mono_ok, poly_ok, PARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wraps run_region_detection with one round of bad-stim replacement.
% A stim is bad in a window if its anchor amp is less negative than
% PARAMS.bad_stim_thresh * (deepest anchor amp across stims), or NaN


%
% SS 2026

reg = run_region_detection(Vb_all, T_VEC, region_chans, stim_idx, ...
    allDepths, t_b_start, t_b_end, mono_win, poly_win, mono_ok, poly_ok, PARAMS);

bad = false(numel(stim_idx), 1);
if reg.mono.peak_detected
    bad = bad | find_bad_slots(reg.mono.stim_peak_amp_anchor, PARAMS.bad_stim_thresh);
end
if reg.poly.peak_detected
    bad = bad | find_bad_slots(reg.poly.stim_peak_amp_anchor, PARAMS.bad_stim_thresh);
end

reg.stim_idx_used = stim_idx(:);
if ~any(bad)
    reg = attach_region_chan_info(reg, region_chans, allDepths);
    return;
end

new_stim_idx = stim_idx(:);
bad_pos      = find(bad);
n_replace    = min(numel(bad_pos), numel(extras_pool));
if n_replace > 0
    new_stim_idx(bad_pos(1:n_replace)) = extras_pool(1:n_replace);
end

reg = run_region_detection(Vb_all, T_VEC, region_chans, new_stim_idx, ...
    allDepths, t_b_start, t_b_end, mono_win, poly_win, mono_ok, poly_ok, PARAMS);

bad2 = false(numel(new_stim_idx), 1);
if reg.mono.peak_detected
    bad2 = bad2 | find_bad_slots(reg.mono.stim_peak_amp_anchor, PARAMS.bad_stim_thresh);
end
if reg.poly.peak_detected
    bad2 = bad2 | find_bad_slots(reg.poly.stim_peak_amp_anchor, PARAMS.bad_stim_thresh);
end

if any(bad2)
    reg.mono = nan_bad_per_stim(reg.mono, bad2);
    reg.poly = nan_bad_per_stim(reg.poly, bad2);
end

reg.stim_idx_used = new_stim_idx(:);
reg = attach_region_chan_info(reg, region_chans, allDepths);

end

% =========================================================================
function reg = attach_region_chan_info(reg, region_chans, allDepths)
% Store the global channel indices and their depths for this region so
% downstream plotting/analysis is self-contained (no need to reload _DATA).
reg.region_chan_global_idx = region_chans(:);
reg.region_chan_depths_um  = double(allDepths(region_chans(:)));
end

% =========================================================================
function bad = find_bad_slots(amps, thresh)
amps    = amps(:);
bad     = isnan(amps) | amps >= 0;
valid   = ~bad;
if ~any(valid); return; end
deepest = min(amps(valid));
bad(valid) = amps(valid) > thresh * deepest;
end

% =========================================================================
function pk = nan_bad_per_stim(pk, bad)
per_stim_fields = { ...
    'mean_amplitude','mean_latency','sem_amplitude', ...
    'mean_amplitude_50pct','mean_latency_50pct', ...
    'mean_amplitude_25pct','mean_latency_25pct', ...
    'n_chans_50pct_per_stim','n_chans_25pct_per_stim', ...
    'stim_peak_amp_anchor','stim_peak_lat_anchor', ...
    'mean_amplitude_3chan','mean_latency_3chan', ...
    'best_channel_amplitude','best_channel_latency', ...
    'all_channel_peaks','all_channel_latencies', ...
    'best_channel_traces', ...
    'trace_per_stim_lowest','trace_per_stim_3chan', ...
    'trace_per_stim_50pct','trace_per_stim_25pct'};
for k = 1:numel(per_stim_fields)
    f = per_stim_fields{k};
    if isfield(pk, f) && ~isempty(pk.(f))
        v = pk.(f);
        v(bad, :) = NaN;
        pk.(f)    = v;
    end
end
end
