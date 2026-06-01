function out = compute_peak_window_times(trace, t_ms, detect_win_ms)
%  Find negative peak and its boundary crossings.
% Locates the most negative point of the input trace and finds boundaries
%  where the trace crosses half-max and 0 (zero crossing) 
% If no crossing is found within the trace, the boundary is set to t_ms(1)
% (rising side) or t_ms(end) (falling side).
%
% SS 2026

out = struct( ...
    'peak_amp_mV',      NaN, ...
    'peak_lat_ms',      NaN, ...
    'halfmax_start_ms', NaN, ...
    'halfmax_end_ms',   NaN, ...
    'zero_start_ms',    NaN, ...
    'zero_end_ms',      NaN, ...
    'valid',            false);

if isempty(trace) || all(isnan(trace)); return; end

%% Peak: minimum within detection window
[~, w_lo] = min(abs(t_ms - detect_win_ms(1)));
[~, w_hi] = min(abs(t_ms - detect_win_ms(2)));
seg = trace(w_lo:w_hi);
if all(isnan(seg)); return; end

[peak_amp, ix_rel] = min(seg);
if peak_amp >= 0; return; end   % no negative peak
peak_idx = w_lo + ix_rel - 1;

out.peak_amp_mV = peak_amp;
out.peak_lat_ms = t_ms(peak_idx);
out.valid       = true;

% Half-max and zero crossings 
out.halfmax_start_ms = find_crossing_left(trace,  t_ms, peak_idx, 0.5 * peak_amp);
out.halfmax_end_ms   = find_crossing_right(trace, t_ms, peak_idx, 0.5 * peak_amp);
out.zero_start_ms    = find_crossing_left(trace,  t_ms, peak_idx, 0);
out.zero_end_ms      = find_crossing_right(trace, t_ms, peak_idx, 0);

end

