function W = override_manual(W, csv_path)
%

if isempty(csv_path) || ~exist(csv_path, 'file')
    error('Manual mode requires manual_window_specifications.csv at: %s', csv_path);
end
M = readtable(csv_path, 'VariableNamingRule', 'preserve');
[rat_M, sess_M, region_M] = norm_keys(M.Rat, M.Session, M.Region);

for i = 1:height(W)
    sel  = strcmp(rat_M, W.Rat{i}) & strcmp(sess_M, W.Sess{i}) & strcmp(region_M, W.Region{i});
    rows = M(sel, :);
    if isempty(rows); continue; end

    if ismember('Mono_Valid', rows.Properties.VariableNames)
        valid = rows(rows.Mono_Valid == 1, :);
        if isempty(valid); valid = rows; end
    else
        valid = rows;
    end
    boundary = max(valid.Mono_End);
    W.Mono_Start(i) = 5;        W.Mono_End(i) = boundary;
    W.Poly_Start(i) = boundary; W.Poly_End(i) = 110;
end
end
