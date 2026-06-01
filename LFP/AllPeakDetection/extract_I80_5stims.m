function [stim_idx, stim_currents_5, extras_pool] = extract_I80_5stims( ...
    I80_DATA, rat_id, sess_name, valid_stims, sorted_current_uA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selects up to 5 valid stims at/above I80 (walk-up); falls back to
% last 5 valid if none above. Also returns extras_pool = additional
% valid stims at/above I80 not in the chosen 5 (ascending), used as
% replacements for bad stims downstream.
%
% SS 2026

stim_idx        = [];
stim_currents_5 = [];
extras_pool     = [];

ratFld = sprintf('Rat%s', rat_id);
if ~isfield(I80_DATA, ratFld);             return; end
if ~isfield(I80_DATA.(ratFld), sess_name); return; end
d = I80_DATA.(ratFld).(sess_name);
if ~isfield(d, 'I80') || isempty(d.I80) || isnan(d.I80); return; end

I80_val      = d.I80;
n_stim_total = numel(sorted_current_uA);
i80_anchor   = find(sorted_current_uA <= I80_val, 1, 'last');
if isempty(i80_anchor); i80_anchor = 1; end

candidate_idx      = i80_anchor:n_stim_total;
valid_above_anchor = candidate_idx(ismember(candidate_idx, valid_stims));

if numel(valid_above_anchor) >= 5
    stim_idx    = valid_above_anchor(1:5);
    extras_pool = valid_above_anchor(6:end);
elseif ~isempty(valid_above_anchor)
    stim_idx = valid_above_anchor;
else
    n_take = min(5, numel(valid_stims));
    if n_take > 0
        vs = valid_stims(:)';
        stim_idx = vs(end - n_take + 1 : end);
    end
end

if isempty(stim_idx); return; end
stim_currents_5 = sorted_current_uA(stim_idx);

end
