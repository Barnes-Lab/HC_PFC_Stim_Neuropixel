function [timeBins, depthBins, allP_norm] = psthByDepthNormalized( ...
    spikeTimes, spikeDepths, global_depthBins, timeBinSize, eventTimes, win, bslWin)
%  Depth-binned PSTH expressed as post-stim minus
% pre-stim firing rate (Hz). For each depth bin and each stim, computes
% baseline rate over bslWin and response rate per time bin over win,
% averages (response - baseline) across stims with any activity. Depth
% bins outside [0.002, 12] Hz mean rate over the whole recording are
% skipped.
%
% INPUT
%   spikeTimes- spike times (s)
%   spikeDepths- per-spike depth (um)
%   global_depthBins - depth bin edges (um)
%   timeBinSize  - bin size (s)
%   eventTimes   -stim times (s)
%   win-response window relative to stim (s)
%   bslWin - baseline window relative to stim (s)
%
% OUTPUT
%   timeBins - bin centers (s)
%   depthBins
%   allP_norm- mean (response - baseline) per depth bin (Hz)
%
% SS 2026 (MODIFIED FROM STEINMETZ)

if isempty(spikeTimes) || isempty(spikeDepths) || isempty(eventTimes)
    timeBins = []; depthBins = []; allP_norm = []; return;
end

min_rate_threshold = 0.002;
max_rate_threshold = 12;

depthBins = global_depthBins;
nD        = length(depthBins) - 1;

edges = win(1):timeBinSize:win(2);
if edges(end) ~= win(2); edges = [edges, win(2)]; end
timeBins = edges(1:end-1) + timeBinSize/2;
nBins    = length(timeBins);

bsl_duration   = diff(bslWin);
total_duration = max(spikeTimes) - min(spikeTimes);
allP_norm      = zeros(nD, nBins);

for d = 1:nD
    inDepth       = spikeDepths > depthBins(d) & spikeDepths <= depthBins(d+1);
    spikesInDepth = spikeTimes(inDepth);
    if isempty(spikesInDepth); continue; end

    depth_rate = length(spikesInDepth) / total_duration;
    if depth_rate < min_rate_threshold || depth_rate > max_rate_threshold; continue; end

    stim_diffs = zeros(length(eventTimes), nBins);
    for e = 1:length(eventTimes)
        baseline_spikes = sum(spikesInDepth >= eventTimes(e) + bslWin(1) & ...
                              spikesInDepth <  eventTimes(e) + bslWin(2));
        baseline_rate   = baseline_spikes / bsl_duration;
        has_activity    = baseline_spikes > 0;

        for t = 1:nBins
            resp_spikes = sum(spikesInDepth >= eventTimes(e) + edges(t) & ...
                              spikesInDepth <  eventTimes(e) + edges(t+1));
            resp_rate   = resp_spikes / timeBinSize;
            if resp_spikes > 0; has_activity = true; end
            stim_diffs(e, t) = resp_rate - baseline_rate;
        end
        if ~has_activity; stim_diffs(e, :) = NaN; end
    end

    valid = stim_diffs(~isnan(stim_diffs(:, 1)), :);
    if ~isempty(valid); allP_norm(d, :) = mean(valid, 1); end
end
end
