function out = NP2_neuron_within_age_posthoc(data, sess_types, metric, recruit_thresh)
% NP2_NEURON_WITHIN_AGE_POSTHOC  6 paired comparisons within each age
% (sessions collapsed):
%   Mono vs Poly (main)   region collapsed
%   PL vs IL (main)       window collapsed
%   Mono: PL vs IL        window = mono fixed
%   Poly: PL vs IL        window = poly fixed
%   PL: Mono vs Poly      region = PL fixed
%   IL: Mono vs Poly      region = IL fixed
%
%
% OUTPUT
%   out- table with Age_Group, Comparison, N_pairs, Mean_Difference, Cohens_d,
%               t_statistic, p_value_raw, p_bonferroni, Significant_p05
%
% SS 2026

ages = {'Y', 'O'};

% pair label, side1 {sess, region, win}, side2 {sess, region, win}
pairs = {
    'Mono vs Poly (main)', {'all','avg','mono'}, {'all','avg','poly'};
    'PL vs IL (main)',     {'all','PL','avg'},   {'all','IL','avg'};
    'Mono: PL vs IL',      {'all','PL','mono'},  {'all','IL','mono'};
    'Poly: PL vs IL',      {'all','PL','poly'},  {'all','IL','poly'};
    'PL: Mono vs Poly',    {'all','PL','mono'},  {'all','PL','poly'};
    'IL: Mono vs Poly',    {'all','IL','mono'},  {'all','IL','poly'};
};

rows = {};
for ai = 1:2
    a = ages{ai};
    rats_list = unique_rats_for_age(data, sess_types, 'PYR', a);
    for p = 1:size(pairs, 1)
        f1 = pairs{p, 2}; f2 = pairs{p, 3};
        v1 = get_per_rat_value(data, sess_types, a, f1{1}, f1{2}, f1{3}, metric, recruit_thresh, rats_list);
        v2 = get_per_rat_value(data, sess_types, a, f2{1}, f2{2}, f2{3}, metric, recruit_thresh, rats_list);

        ok = ~isnan(v1) & ~isnan(v2);
        vA = v1(ok); vB = v2(ok);
        n  = numel(vA);
        if n < 2
            rows(end+1, :) = {a, pairs{p, 1}, n, NaN, NaN, NaN, NaN};
            continue;
        end
        [~, pval, ~, st] = ttest(vA, vB);
        rows(end+1, :) = {a, pairs{p, 1}, n, mean(vA - vB), ...
            NP2_cohens_d(vA, vB, 'paired'), st.tstat, pval};
    end
end

n_comp = size(pairs, 1);
out = cell2table(rows, 'VariableNames', {'Age_Group','Comparison','N_pairs', ...
    'Mean_Difference','Cohens_d','t_statistic','p_value_raw'});
out.p_bonferroni    = min(out.p_value_raw * n_comp, 1);
out.Significant_p05 = out.p_value_raw < 0.05;
end
