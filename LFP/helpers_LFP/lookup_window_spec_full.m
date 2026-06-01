function [mono_win, poly_win, mono_valid, poly_valid, has_override] = ...
    lookup_window_spec_full(windowTbl, rat_id, sess_name, shank, region, ...
                            default_mono_win, default_poly_win)
% Look up per-rat per-session per-shank per-region MONO and POLY windows
% from manual_window_specifications.csv. 
%
% CSV format
%   Rat, Session, Shank, Region, Mono_Start, Mono_End, Poly_Start,
%   Poly_End, Mono_Valid, Poly_Valid, Mono_BestChan, Poly_BestChan
%
% SS 2026

mono_win     = default_mono_win;
poly_win     = default_poly_win;
mono_valid   = true;
poly_valid   = true;
has_override = false;

if isempty(windowTbl); return; end

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
poly_win     = [row.Poly_Start, row.Poly_End];
mono_valid   = logical(row.Mono_Valid);
poly_valid   = logical(row.Poly_Valid);
has_override = true;

end
