function W = override_fepsp(W, csv_path, win_type)
% OVERRIDE_FEPSP  
% SS 2026

if isempty(csv_path) || ~exist(csv_path, 'file')
    error('fepsp_peak mode requires peak_window_times_post_i80.csv at: %s', csv_path);
end
F = readtable(csv_path, 'VariableNamingRule', 'preserve');
[rat_F, sess_F, region_F] = norm_keys(F.Rat, F.Session, F.Region);
win_F = strtrim(string(F.Window));

switch lower(win_type)
    case 'tight';      cs = 'HalfMax_Start_ms'; ce = 'HalfMax_End_ms';
    case 'window_all'; cs = 'Zero_Start_ms';    ce = 'Zero_End_ms';
    otherwise; error('Unknown fepsp_win_type: %s (use ''tight'' or ''window_all'')', win_type);
end

for i = 1:height(W)
    base  = strcmp(rat_F, W.Rat{i}) & strcmp(sess_F, W.Sess{i}) & strcmp(region_F, W.Region{i});
    sel_m = base & strcmpi(win_F, 'Mono');
    sel_p = base & strcmpi(win_F, 'Poly');
    if any(sel_m)
        W.Mono_Start(i) = F.(cs)(sel_m);
        W.Mono_End(i)   = F.(ce)(sel_m);
    end
    if any(sel_p)
        W.Poly_Start(i) = F.(cs)(sel_p);
        W.Poly_End(i)   = F.(ce)(sel_p);
    end
end
end
