function [event_times, currents_5, method] = select_I80_event_times_neuron( ...
    I80_DATA, rat_id, sess_name, T_sess, cur_col)
%Spike-data version  of extract_I80_5stims.
% Selects 5 stim event times based on chosen-shank I80 from I80_DATA.
%

% INPUT
%   I80_DATAloaded from I80_DATA.mat (LFP generated)
%   rat_id ; sess_name cur_col-stim current col    (ALL CHAR)
%   T_sess -ALL_STIM_TIMES 
%
% OUTPUT
%   event_times -Time_s for selected stims 
%   currents_5  uA value for stims

%
% SS 2026

event_times = []; currents_5 = []; method = 'no_data';
if isempty(T_sess); return; end


keep_col = '';
if     ismember('Keep', T_sess.Properties.VariableNames); keep_col = 'Keep';
elseif ismember('keep', T_sess.Properties.VariableNames); keep_col = 'keep';
end

T = sortrows(T_sess, cur_col);
n = height(T);
if ~isempty(keep_col)
    is_valid = ~strcmpi(string(T.(keep_col)), 'N');
else
    is_valid = true(n, 1);
end

% I80 lookup
ratFld  = sprintf('Rat%s', rat_id);
I80_val = NaN;
if isfield(I80_DATA, ratFld) && isfield(I80_DATA.(ratFld), sess_name)
    d = I80_DATA.(ratFld).(sess_name);
    if isfield(d, 'I80') && ~isempty(d.I80) && ~isnan(d.I80)
        I80_val = d.I80;
    end
end

currs = T.(cur_col);
sel   = [];
if ~isnan(I80_val)
    anchor = find(currs <= I80_val, 1, 'last');
    if isempty(anchor); anchor = 1; end
    candidate_pos = anchor : n;
    valid_above   = candidate_pos(is_valid(candidate_pos));
    if numel(valid_above) >= 5
        sel = valid_above(1:5); method = 'I80';
    elseif ~isempty(valid_above)
        sel = valid_above;       method = 'I80_partial';
    else
        method = 'last_valid';
    end
else
    method = 'no_I80';
end

% Fallback: last 5 valid
if isempty(sel)
    valid_all = find(is_valid);
    if isempty(valid_all); return; end
    n_take = min(5, numel(valid_all));
    sel = valid_all(end - n_take + 1 : end);
    if ~strcmp(method, 'I80_partial'); method = 'last_valid'; end
end

event_times = T.Time_s(sel);
currents_5  = T.(cur_col)(sel);
end
