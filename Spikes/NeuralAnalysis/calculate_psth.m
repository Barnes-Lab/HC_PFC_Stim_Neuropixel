function psth = calculate_psth(spike_times, event_times, time_bins)
% CALCULATE_PSTH  Per-trial-averaged spike rate (Hz) per time bin, relative
% to event times.
%
% INPUT
%   spike_times(s)
%   event_times(s)
%   time_bins bin edges  (s)
%
% OUTPUT
%   psth  -mean firing rate per bin (Hz)
%
% Adapted from Steinmetz github repo

n_bins = length(time_bins) - 1;
psth   = zeros(1, n_bins);

for e = 1:length(event_times)
    rel_spikes = spike_times - event_times(e);
    for b = 1:n_bins
        psth(b) = psth(b) + sum(rel_spikes >= time_bins(b) & rel_spikes < time_bins(b+1));
    end
end

psth = psth / length(event_times) / diff(time_bins(1:2));
end
