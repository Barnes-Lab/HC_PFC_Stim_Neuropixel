function [mono_win, mono_valid, has_override] = lookup_window_spec( ...
    windowTbl, rat_id, sess_name, shank, region, default_mono_win)
%
% Look up per-rat per-session per-shank per-region mono window from
% manual_window_specifications.csv. Returns default + has_override=false
% if no row matches.

% SS 2026

mono_win     = default_mono_win;
mono_valid   = true;
has_override = false;

if isempty(windowTbl); return; end

%% Trim string columns (CSV had spaces like 'iHC    ')
rat_col  = strtrim(string(windowTbl.Rat));
sess_col = strtrim(string(windowTbl.Session));
reg_col  = strtrim(string(windowTbl.Region));

mask = (rat_col  == string(rat_id))    & ...
       (sess_col == string(sess_name)) & ...
       (windowTbl.Shank == shank)      & ...
       (reg_col  == string(region));

if ~any(mask); return; end

row = windowTbl(find(mask, 1), :);

mono_win     = [row.Mono_Start, row.Mono_End];
mono_valid   = logical(row.Mono_Valid);
has_override = true;

end
