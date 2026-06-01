function [is_walled, frac_walled] = check_window_wall_hit(pk, mono_window_ms, ...
    wall_tol_ms, wall_hit_frac)
% Looks for cases where peak is detected at window edge for many stims
% indicates incorrect detection
%
% A stim is "walled" if its per-stim peak latency falls within wall_tol_ms
% of mono_window_ms(2). is_walled = (frac walled) >= wall_hit_frac.

%
% SS 2026

is_walled   = false;
frac_walled = NaN;

if isempty(pk) || ~isstruct(pk) || ~isfield(pk, 'peak_detected') || ~pk.peak_detected
    return;
end

lats  = pk.best_channel_latency;
valid = ~isnan(lats);
if sum(valid) == 0; return; end

frac_walled = sum(lats(valid) >= (mono_window_ms(2) - wall_tol_ms)) / sum(valid);
is_walled   = frac_walled >= wall_hit_frac;

end
