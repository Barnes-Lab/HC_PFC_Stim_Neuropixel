function pk = detect_mono_v2_findpeaks(V_proc, t, channels, valid_stims, ...
    candidate_local, allDepths, mono_window_ms, t_b_start, t_b_end, varargin)
%Detects mono peak using a template
% Algorithm:
%   1. Per-channel template = mean across N strongest stims.
%   2. Per-channel sigma from baseline window.
%   3. findpeaks on -template per candidate channel within mono window:
%        MinPeakProminence = PRM_prominence_sigma * sigma_c
%        MinPeakDistance, MinPeakWidth
%   4. For each candidate channel, take the MOST negative qualifying peak.
%   5. Best channel = candidate with the most negative qualifying peak.
%   6. Per-stim extraction in [template_lat +/- PRM_per_trial_half_win_ms]
%      (template-anchored, NOT walled by mono_window).
%
% SS 2026

%% PARAMETERS TO CHANGE
PRM_n_strongest           = 3;
PRM_prominence_sigma      = 3;
PRM_min_peak_distance_ms  = 5;
PRM_min_peak_width_ms     = 0.5;
PRM_per_trial_half_win_ms = 10;
Extract_varargin;

n_stim_v = size(V_proc, 1);
n_chan_r = size(V_proc, 2);
n_t      = size(V_proc, 3);

pk = empty_peak_struct(n_stim_v, n_chan_r, mono_window_ms);

if numel(candidate_local) < 3 || n_stim_v < PRM_n_strongest; return; end

%% Strong stims channels
[~, mw_lo] = min(abs(t - mono_window_ms(1)));
[~, mw_hi] = min(abs(t - mono_window_ms(2)));
per_stim_min = nan(n_stim_v, 1);
for i = 1:n_stim_v
    seg = squeeze(V_proc(i, candidate_local, mw_lo:mw_hi));
    per_stim_min(i) = min(seg(:));
end
[~, order]       = sort(per_stim_min, 'ascend');
strong_idx_local = order(1:PRM_n_strongest);
strong_stims     = valid_stims(strong_idx_local);

%% Templates per channel
templates = squeeze(mean(V_proc(strong_idx_local, :, :), 1));   % [nChanReg x nT]
sigmas    = std(templates(:, t_b_start:t_b_end), 0, 2);

%% findpeaks per candidate
search_idx     = mw_lo:mw_hi;
t_search       = t(search_idx);
dt_ms          = t(2) - t(1);
min_dist_samp  = max(1, round(PRM_min_peak_distance_ms / dt_ms));
min_width_samp = max(1, round(PRM_min_peak_width_ms    / dt_ms));

cand_amp = nan(n_chan_r, 1);
cand_lat = nan(n_chan_r, 1);

for ci = candidate_local(:)'   % LOOP: per-channel findpeaks
    if sigmas(ci) <= 0 || ~isfinite(sigmas(ci)); continue; end
    [pks_pos, locs] = findpeaks(-templates(ci, search_idx), ...
        'MinPeakProminence', PRM_prominence_sigma * sigmas(ci), ...
        'MinPeakDistance',   min_dist_samp, ...
        'MinPeakWidth',      min_width_samp);
    if isempty(pks_pos); continue; end
    [a, k]       = min(-pks_pos(:));
    cand_amp(ci) = a;
    cand_lat(ci) = t_search(locs(k));
end

if all(isnan(cand_amp)); return; end

%% Best candidate channel
[anchor_amp, anchor_local] = min(cand_amp);
if isnan(anchor_amp) || anchor_amp >= 0; return; end
anchor_lat    = cand_lat(anchor_local);
anchor_global = channels(anchor_local);

%% Anatomical neighbours
if anchor_local == 1
    neighbours_local = [2; 3];
elseif anchor_local == n_chan_r
    neighbours_local = [n_chan_r - 1; n_chan_r - 2];
else
    neighbours_local = [anchor_local - 1; anchor_local + 1];
end
neighbours_local = neighbours_local(neighbours_local >= 1 & neighbours_local <= n_chan_r);
if numel(neighbours_local) < 2; return; end
three_local  = [anchor_local; neighbours_local(:)];
three_global = channels(three_local);

%% Per-stim extraction 
[~, pw_lo] = min(abs(t - (anchor_lat - PRM_per_trial_half_win_ms)));
[~, pw_hi] = min(abs(t - (anchor_lat + PRM_per_trial_half_win_ms)));
pw_lo = max(pw_lo, 1); pw_hi = min(pw_hi, n_t);
t_pw  = t(pw_lo:pw_hi);

all_amp     = nan(n_stim_v, n_chan_r);
all_lat     = nan(n_stim_v, n_chan_r);
best_traces = nan(n_stim_v, n_t);
for i = 1:n_stim_v
    V_seg          = squeeze(V_proc(i, :, pw_lo:pw_hi));
    [a, ix]        = min(V_seg, [], 2);
    all_amp(i, :)  = a';
    all_lat(i, :)  = t_pw(ix);
    best_traces(i, :) = V_proc(i, anchor_local, :);
end
all_amp(all_amp >= 0) = 0;

best_amp_per_trial = all_amp(:, anchor_local);
best_lat_per_trial = all_lat(:, anchor_local);
mean_amp = mean(all_amp(:, three_local), 2, 'omitnan');
mean_lat = mean(all_lat(:, three_local), 2, 'omitnan');
sem_amp  = std(all_amp(:, three_local), 0, 2, 'omitnan') ./ sqrt(numel(three_local));

stim_rank = (1:n_stim_v)';
keep      = ~isnan(best_amp_per_trial);
if sum(keep) >= 4
    [rho, pval] = corr(stim_rank(keep), best_amp_per_trial(keep), 'type', 'Spearman');
else
    rho = NaN; pval = NaN;
end

pk.peak_detected            = true;
pk.detection_method         = 'mono_v2_findpeaks';
pk.detection_signal         = 'findpeaks_per_channel';
pk.best_channel             = anchor_global;
pk.best_channel_depth       = allDepths(anchor_global);
pk.template_latency         = anchor_lat;
pk.template_amplitude       = anchor_amp;
pk.search_window            = mono_window_ms;
pk.strong_stims             = strong_stims(:);
pk.threshold_channels       = three_global;
pk.n_threshold_channels     = numel(three_global);
pk.best_channel_amplitude   = best_amp_per_trial;
pk.best_channel_latency     = best_lat_per_trial;
pk.mean_amplitude           = mean_amp;
pk.mean_latency             = mean_lat;
pk.sem_amplitude            = sem_amp;
pk.all_channel_peaks        = all_amp;
pk.all_channel_latencies    = all_lat;
pk.best_channel_traces      = best_traces;
pk.monotonicity_rho         = rho;
pk.monotonicity_p           = pval;
pk.excluded_channels        = setdiff(channels, channels(candidate_local));

end
% Change this each time so setting this inline
function p = empty_peak_struct(n_stim, n_chan, lat_range)
p = struct();
p.peak_detected             = false;
p.detection_method          = 'none';
p.detection_signal          = 'findpeaks_per_channel';
p.best_channel              = NaN;
p.best_channel_depth        = NaN;
p.template_latency          = NaN;
p.template_amplitude        = NaN;
p.search_window             = lat_range;
p.strong_stims              = [];
p.threshold_channels        = [];
p.n_threshold_channels      = 0;
p.best_channel_amplitude    = nan(n_stim, 1);
p.best_channel_latency      = nan(n_stim, 1);
p.mean_amplitude            = nan(n_stim, 1);
p.mean_latency              = nan(n_stim, 1);
p.sem_amplitude             = nan(n_stim, 1);
p.all_channel_peaks         = nan(n_stim, n_chan);
p.all_channel_latencies     = nan(n_stim, n_chan);
p.best_channel_traces       = [];
p.monotonicity_rho          = NaN;
p.monotonicity_p            = NaN;
p.excluded_channels         = [];
end
