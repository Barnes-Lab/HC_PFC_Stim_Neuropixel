function out = NP2_neuron_between_age_by_session_posthoc(data, sess_types, metric, recruit_thresh)
% NP2_NEURON_BETWEEN_AGE_BY_SESSION_POSTHOC  Y vs O within each session.
% 7 conditions per session: PL-mono, PL-poly, IL-mono, IL-poly,
% Mono-all (across channels), Poly-all All-All
% and windows pooled). Total = 7 x n_sess rows.
% Correected with bonferoni 
%%


comps = {
    'mono-PL',   'PL',  'mono';
    'poly-PL',   'PL',  'poly';
    'mono-IL',   'IL',  'mono';
    'poly-IL',   'IL',  'poly';
    'mono-all',  'avg', 'mono';
    'poly-all',  'avg', 'poly';
    'all-all',   'avg', 'avg';
};

rats_Y = unique_rats_for_age(data, sess_types, 'PYR', 'Y');
rats_O = unique_rats_for_age(data, sess_types, 'PYR', 'O');

n_per = size(comps, 1);
rows  = cell(numel(sess_types) * n_per, 11);
k = 0;
for s = 1:numel(sess_types)
    sess = sess_types{s};
    for ci = 1:n_per
        cond   = comps{ci, 1};
        region = comps{ci, 2};
        win    = comps{ci, 3};
        yY = get_per_rat_value(data, sess_types, 'Y', sess, region, win, metric, recruit_thresh, rats_Y);
        yO = get_per_rat_value(data, sess_types, 'O', sess, region, win, metric, recruit_thresh, rats_O);
        yY = yY(~isnan(yY)); yO = yO(~isnan(yO));

        k = k + 1;
        if numel(yY) < 2 || numel(yO) < 2
            rows(k, :) = {sess, cond, numel(yY), numel(yO), NaN, NaN, NaN, NaN, NaN, NaN, NaN};
            continue;
        end
        [~, pval, ~, st] = ttest2(yY, yO);
        rows(k, :) = {sess, cond, numel(yY), numel(yO), mean(yY), mean(yO), ...
            mean(yY) - mean(yO), NP2_cohens_d(yY, yO, 'unpaired'), ...
            st.tstat, st.df, pval};
    end
end

out = cell2table(rows, 'VariableNames', {'Sess','Condition','N_Young','N_Old', ...
    'Young_Mean','Old_Mean','Mean_Difference','Cohens_d','t_statistic','df','p_value_raw'});
out.p_bonferroni    = min(out.p_value_raw * n_per, 1);
out.Significant_p05 = out.p_value_raw < 0.05;
end
