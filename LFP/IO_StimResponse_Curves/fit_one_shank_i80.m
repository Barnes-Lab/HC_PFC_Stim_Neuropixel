function sr = fit_one_shank_i80(shank_struct, shank_idx, sess_idx, ...
                                 stim_tbl, rat_id, varargin)
% Pulls IL mono amplitudes for one (shank, session) from MONO_DATA,
% and runs cascasde to fit curves

%
% Example call:
%   sr = fit_one_shank_i80(MONO_DATA.Rat10987.Shank3, 3, 1, ...
%        stim_tbl, '10987', 'PRM_force_method', 'exp_sat');
%
% 
% DEPENDS  try_fit_cascade.m
%
% SS 2026

%% PARAMETERS
PRM_R2_min         = 0.50;
PRM_R2_forced      = 0.20;
PRM_I80_range      = [50, 900];
PRM_flat_thresh_mV = 0.10;
PRM_amp_field      = 'mean_amplitude';
PRM_force_method   = '';
Extract_varargin;

sr = empty_shank_result();
sr.shank_idx = shank_idx;

if isfield(shank_struct, 'layer_designation')
    sr.layer = shank_struct.layer_designation;
end

if ~isfield(shank_struct, 'sessions') || numel(shank_struct.sessions) < sess_idx
    return;
end

sd = shank_struct.sessions(sess_idx);
sr.sess_name         = sd.name;
sr.stim_currents_all = sd.stim_currents(:);
sr.valid_stims       = sd.valid_stims(:);

if ~isstruct(sd.IL.mono) || ~sd.IL.mono.peak_detected
    sr.Status = 'no_peak'; return;
end

if isfield(sd.IL, 'channel_mode')
    sr.channel_mode = sd.IL.channel_mode;
end

if ~isfield(sd.IL.mono, PRM_amp_field)
    sr.Status = 'amp_field_missing'; return;
end

amps  = sd.IL.mono.(PRM_amp_field); amps = amps(:);
x_all = sr.stim_currents_all(sr.valid_stims);

if numel(x_all) ~= numel(amps)
    sr.Status = 'amp_len_mismatch'; return;
end

keep       = ~isnan(x_all) & ~isnan(amps) & isfinite(x_all) & isfinite(amps);
sr.x       = x_all(keep);
sr.y       = amps(keep);
valid_kept = sr.valid_stims(keep);

if numel(sr.x) < 4
    sr.Status = 'n_lt_4'; return;
end

sr.has_data = true;
sr.is_flat  = (max(sr.y) - min(sr.y)) < PRM_flat_thresh_mV;

maybe_mask = build_maybe_mask_silent(stim_tbl, rat_id, sr.sess_name, valid_kept);

sr.fit = try_fit_cascade(sr.x, sr.y, maybe_mask, ...
    'PRM_R2_min',       PRM_R2_min, ...
    'PRM_R2_forced',    PRM_R2_forced, ...
    'PRM_I80_range',    PRM_I80_range, ...
    'PRM_force_method', PRM_force_method);

sr.Status = resolve_status(sr.fit, PRM_R2_min, sr.is_flat);

end

% Inline
function maybe_mask = build_maybe_mask_silent(stim_tbl, rat_id, sess_name, valid_kept)
% Same as build_maybe_mask.m but no console output (called per-shank).
maybe_mask = false(numel(valid_kept), 1);
if isempty(stim_tbl); return; end
rows = stim_tbl( ...
    strcmp(string(stim_tbl.Rat),     string(rat_id)) & ...
    strcmp(string(stim_tbl.Session), string(sess_name)), :);
if isempty(rows); return; end
if ~any(strcmp(rows.Properties.VariableNames, 'Maybe')); return; end
maybe_col = strip(string(rows.Maybe));
for j = 1:numel(valid_kept)
    idx = valid_kept(j);
    if idx <= numel(maybe_col) && strcmpi(maybe_col(idx), 'N')
        maybe_mask(j) = true;
    end
end
end

% Inlines
function status = resolve_status(fit, R2_thr, is_flat)
m = fit.Method;
if strcmp(m, 'none')
    status = 'failed';
elseif strcmp(m, 'data_driven')
    if is_flat; status = 'flat_data_driven'; else; status = 'data_driven'; end
elseif endsWith(m, 'forced_method')
    status = 'ok_force_method';
elseif endsWith(m, '_forced_no_maybe')
    status = 'ok_forced_no_maybe';
elseif endsWith(m, '_forced')
    status = 'ok_forced';
elseif endsWith(m, '_no_maybe')
    status = 'ok_no_maybe';
elseif ~isnan(fit.R2) && fit.R2 >= R2_thr
    status = 'ok';
else
    status = 'below_R2';
end
end
