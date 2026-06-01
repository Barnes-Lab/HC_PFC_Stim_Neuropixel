function out = NP2_neuron_between_age_posthoc(data, sess_types, metric, recruit_thresh)
% NP2_NEURON_BETWEEN_AGE_POSTHOC  4 between-age comparisons (Y vs O,
% sessions collapsed): PL-mono, PL-poly, IL-mono, IL-poly.
%
% Mathematical rationale
%   Two-sample t (pooled SD) on rat-level values; df = n_Y + n_O - 2.
%   Cohen's d with pooled SD. Bonferroni across the 4 cells.
%
% INPUT
%   data            struct
%   sess_types      cell
%   metric          char    (see get_per_rat_value)
%   recruit_thresh  scalar  Hz
%
% OUTPUT
%   out  table  Comparison, N_Young, N_Old, Young_Mean, Old_Mean,
%               Mean_Difference, Cohens_d, t_statistic, df, p_value_raw,
%               p_bonferroni, Significant_p05
%
% SS 2026

regions = {'PL', 'IL'};
windows = {'mono', 'poly'};

rats_Y = unique_rats_for_age(data, sess_types, 'PYR', 'Y');
rats_O = unique_rats_for_age(data, sess_types, 'PYR', 'O');

rows = {};
for ri = 1:2
    for wi = 1:2
        cond = sprintf('%s-%s', regions{ri}, windows{wi});
        yY = get_per_rat_value(data, sess_types, 'Y', 'all', regions{ri}, windows{wi}, metric, recruit_thresh, rats_Y);
        yO = get_per_rat_value(data, sess_types, 'O', 'all', regions{ri}, windows{wi}, metric, recruit_thresh, rats_O);
        yY = yY(~isnan(yY)); yO = yO(~isnan(yO));

        if numel(yY) < 2 || numel(yO) < 2
            rows(end+1, :) = {cond, numel(yY), numel(yO), NaN, NaN, NaN, NaN, NaN, NaN, NaN};
            continue;
        end
        [~, pval, ~, st] = ttest2(yY, yO);
        rows(end+1, :) = {cond, numel(yY), numel(yO), mean(yY), mean(yO), ...
            mean(yY) - mean(yO), NP2_cohens_d(yY, yO, 'unpaired'), ...
            st.tstat, st.df, pval};
    end
end

n_comp = 4;
out = cell2table(rows, 'VariableNames', {'Comparison','N_Young','N_Old', ...
    'Young_Mean','Old_Mean','Mean_Difference','Cohens_d','t_statistic','df','p_value_raw'});
out.p_bonferroni    = min(out.p_value_raw * n_comp, 1);
out.Significant_p05 = out.p_value_raw < 0.05;
end
