function [out, det_win] = compute_for_layer(rat_entry, shIdx, sess_name, ...
                                            region, win_fld, T_VEC)
% Pull lowest-channel trace per layer
% Compute peak/crossings 
out = struct('peak_amp_mV', NaN, 'peak_lat_ms', NaN, ...
    'halfmax_start_ms', NaN, 'halfmax_end_ms', NaN, ...
    'zero_start_ms',    NaN, 'zero_end_ms',    NaN, 'valid', false);
det_win = [NaN NaN];

if isnan(shIdx); return; end
shFld = sprintf('Shank%d', shIdx);
if ~isfield(rat_entry, shFld); return; end
if ~isfield(rat_entry.(shFld), 'sessions'); return; end

sess_struct = rat_entry.(shFld).sessions;
sess_idx = NaN;
for s = 1:numel(sess_struct)
    if isfield(sess_struct(s), 'name') && strcmp(sess_struct(s).name, sess_name)
        sess_idx = s; break;
    end
end
if isnan(sess_idx); return; end

if ~isfield(sess_struct(sess_idx), region); return; end
reg = sess_struct(sess_idx).(region);
if ~isfield(reg, win_fld); return; end
pk = reg.(win_fld);
if isempty(pk) || ~pk.peak_detected; return; end
if ~isfield(pk, 'trace_per_stim_lowest') || isempty(pk.trace_per_stim_lowest)
    return;
end

ts = pk.trace_per_stim_lowest;             % [K x nT]
ok = ~all(isnan(ts), 2);
if ~any(ok); return; end
trace = mean(ts(ok, :), 1, 'omitnan');     % [1 x nT]

det_win = pk.search_window;
out = compute_peak_window_times(trace, T_VEC, det_win);
end


