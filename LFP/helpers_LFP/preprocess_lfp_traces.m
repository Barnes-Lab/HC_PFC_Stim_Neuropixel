function out = preprocess_lfp_traces(Vb_all, t, channels, valid_stims, ...
    t_b_start, t_b_end, varargin)
% Region preprocessing steps:
%   1) Pre-stim baseline subtraction at [t_b_start, t_b_end]
%   2) Pre-peak baseline subtraction at PRM_prepeak_baseline_ms
%   3) Savitzky-Golay smoothing (order 3, PRM_smooth_window_ms)
%   4) Drift exclusion: channels whose mean pre-peak value across the
%      PRM_n_strongest  exceeds PRM_drift_thresh_mV are removed
%
%
% SS 2026

%% PARAMETERS TO CHANGE
PRM_prepeak_baseline_ms = [5, 9];
PRM_smooth_window_ms    = 3;
PRM_drift_thresh_mV     = 0.15;
PRM_strong_window_ms    = [10, 40];
PRM_n_strongest         = 3;
Extract_varargin;
%%
n_stim_v = numel(valid_stims);
n_chan_r = numel(channels);
n_t      = numel(t);

[~, pp_lo] = min(abs(t - PRM_prepeak_baseline_ms(1)));
[~, pp_hi] = min(abs(t - PRM_prepeak_baseline_ms(2)));
[~, sw_lo] = min(abs(t - PRM_strong_window_ms(1)));
[~, sw_hi] = min(abs(t - PRM_strong_window_ms(2)));

dt        = t(2) - t(1);
sm_n      = round(PRM_smooth_window_ms / dt);
if mod(sm_n, 2) == 0; sm_n = sm_n + 1; end
sm_n      = max(5, sm_n);

V_proc        = zeros(n_stim_v, n_chan_r, n_t);
prepeak_means = zeros(n_stim_v, n_chan_r);

for i = 1:n_stim_v   % -per-stim baseline correction differs across stims
    V_stim = squeeze(Vb_all(valid_stims(i), :, :));
    V_reg  = V_stim(channels, :);
    V_reg  = V_reg - mean(V_reg(:, t_b_start:t_b_end), 2);
    prepeak_means(i, :) = mean(V_reg(:, pp_lo:pp_hi), 2);
    V_reg  = V_reg - prepeak_means(i, :)';
    V_reg  = sgolayfilt(V_reg, 3, sm_n, [], 2);
    V_proc(i, :, :) = V_reg;
end

%% Strong stims ie N most-negative across-channel mins 
per_stim_min = nan(n_stim_v, 1);
for i = 1:n_stim_v
    seg = squeeze(V_proc(i, :, sw_lo:sw_hi));
    per_stim_min(i) = min(seg(:));
end
n_pick           = min(PRM_n_strongest, n_stim_v);
[~, order]       = sort(per_stim_min, 'ascend');
strong_idx_local = order(1:n_pick);

%% Drift exclusion
mean_prepeak_strong = mean(prepeak_means(strong_idx_local, :), 1);
drifting_local      = abs(mean_prepeak_strong) > PRM_drift_thresh_mV;
candidate_local     = find(~drifting_local);

out  = struct();
out.V_proc   = V_proc;
out.candidate_local  = candidate_local;
out.excluded_global  = channels(drifting_local);
out.strong_idx_local = strong_idx_local;
out.strong_stims  = valid_stims(strong_idx_local);

end
