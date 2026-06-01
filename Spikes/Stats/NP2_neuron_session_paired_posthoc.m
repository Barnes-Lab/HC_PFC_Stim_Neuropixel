function out = NP2_neuron_session_paired_posthoc(data, sess_types, metric, recruit_thresh)
% NP2_NEURON_SESSION_PAIRED_POSTHOC  iHC vs vHC paired t-test per rat.
% Per-rat session value = rat's mean over 
%
% SS 2026

rats_Y = unique_rats_for_age(data, sess_types, 'PYR', 'Y');
rats_O = unique_rats_for_age(data, sess_types, 'PYR', 'O');

vY_iHC = get_per_rat_value(data, sess_types, 'Y', 'iHC', 'avg', 'avg', metric, recruit_thresh, rats_Y);
vY_vHC = get_per_rat_value(data, sess_types, 'Y', 'vHC', 'avg', 'avg', metric, recruit_thresh, rats_Y);
vO_iHC = get_per_rat_value(data, sess_types, 'O', 'iHC', 'avg', 'avg', metric, recruit_thresh, rats_O);
vO_vHC = get_per_rat_value(data, sess_types, 'O', 'vHC', 'avg', 'avg', metric, recruit_thresh, rats_O);

groups = {'All', 'Y', 'O'};
rows   = {};
for g = 1:3
    switch groups{g}
        case 'All'; vi = [vY_iHC; vO_iHC]; vv = [vY_vHC; vO_vHC];
        case 'Y';   vi = vY_iHC;            vv = vY_vHC;
        case 'O';   vi = vO_iHC;            vv = vO_vHC;
    end
    ok = ~isnan(vi) & ~isnan(vv);
    vi = vi(ok); vv = vv(ok);
    n  = numel(vi);
    if n < 2
        rows(end+1, :) = {groups{g}, n, NaN, NaN, NaN, NaN};
        continue;
    end
    [~, pval, ~, st] = ttest(vi, vv);
    rows(end+1, :) = {groups{g}, n, mean(vi - vv), ...
        NP2_cohens_d(vi, vv, 'paired'), st.tstat, pval};
end

n_comp = 3;
out = cell2table(rows, 'VariableNames', {'Group','N_pairs', ...
    'Mean_Difference','Cohens_d','t_statistic','p_value_raw'});
out.p_bonferroni    = min(out.p_value_raw * n_comp, 1);
out.Significant_p05 = out.p_value_raw < 0.05;
end
