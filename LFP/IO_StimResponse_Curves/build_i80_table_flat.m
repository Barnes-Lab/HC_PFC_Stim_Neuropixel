function T = build_i80_table_flat(I80_DATA, RATS, AGES)
%
% Flattens I80_DATA (nested struct, all 4 shanks per rat x session)
% into a single-row-per-(rat, session) table containing the chosen
% shank's fit
%
% cols
%   Rat, Age, Stim_Site, Layer, Shank, Session,
%   I80, Amin, Amax, I50, k, R2, Method, Status,
%   Chosen_Source, L5_I80, Best_R2_Shank, L5_vs_Best_diff_uA
%
%
% SS 2026

template = struct( ...
    'Rat','', 'Age','', 'Stim_Site','', 'Layer','', ...
    'Shank',NaN, 'Session',NaN, ...
    'I80',NaN, 'Amin',NaN, 'Amax',NaN, 'I50',NaN, 'k',NaN, 'R2',NaN, ...
    'Method','', 'Status','', 'Chosen_Source','', ...
    'L5_I80',NaN, 'Best_R2_Shank',NaN, 'L5_vs_Best_diff_uA',NaN);

rows = repmat(template, 2 * numel(RATS), 1);
irow = 0;

for iRat = 1:numel(RATS)
    rat_id = RATS{iRat};
    age    = AGES{iRat};
    ratFld = sprintf('Rat%s', rat_id);
    if ~isfield(I80_DATA, ratFld); continue; end

    sess_names = fieldnames(I80_DATA.(ratFld));
    for s = 1:numel(sess_names)
        sn = sess_names{s};
        d  = I80_DATA.(ratFld).(sn);
        chosen = d.chosen_shank;
        if isnan(chosen) || chosen < 1 || chosen > numel(d.shank_results); continue; end

        sr = d.shank_results(chosen);
        sess_idx = find(strcmp({'iHC','vHC'}, sn), 1);

        irow = irow + 1;
        rows(irow).Rat                = rat_id;
        rows(irow).Age                = age;
        rows(irow).Stim_Site          = sn;
        rows(irow).Layer              = sr.layer;
        rows(irow).Shank              = chosen;
        rows(irow).Session            = sess_idx;
        rows(irow).I80                = d.I80;
        rows(irow).Amin               = sr.fit.Amin;
        rows(irow).Amax               = sr.fit.Amax;
        rows(irow).I50                = sr.fit.I50;
        rows(irow).k                  = sr.fit.k;
        rows(irow).R2                 = sr.fit.R2;
        rows(irow).Method             = sr.fit.Method;
        rows(irow).Status             = sr.Status;
        rows(irow).Chosen_Source      = d.chosen_source;
        rows(irow).L5_I80             = d.L5_I80;
        rows(irow).Best_R2_Shank      = d.Best_R2_Shank;
        rows(irow).L5_vs_Best_diff_uA = d.L5_vs_Best_diff_uA;
    end
end

T = struct2table(rows(1:irow));

end
