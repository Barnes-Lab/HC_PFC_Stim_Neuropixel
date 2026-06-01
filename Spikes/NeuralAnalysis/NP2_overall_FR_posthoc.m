function out = NP2_overall_FR_posthoc(data, sess_types, mode, recruit_thresh)
% Y vs O two-sample t on per-rat win_combined dFR
% (PYR), for each (region, sess_filter) 
% Modes:
%   'all'+per-rat mean - all PYR neurons
%   'responders'  per-rat mean -responders (per win > threshold)
% USED IN REPORTING ONLY
% SS 2026

regions       = {'PL', 'IL', 'all'};
region_titles = {'PL', 'IL', 'PL+IL'};
sess_filters  = {'all', 'iHC', 'vHC'};
sess_titles   = {'Combined', 'iHC', 'vHC'};

resp_only = strcmp(mode, 'responders');
rats_Y = unique_rats_for_age(data, sess_types, 'PYR', 'Y');
rats_O = unique_rats_for_age(data, sess_types, 'PYR', 'O');

rows = {};
for ri = 1:3
    for sj = 1:3
        vY = aggregate_win_combined_per_rat(data, sess_types, 'Y', ...
            sess_filters{sj}, regions{ri}, resp_only, recruit_thresh, rats_Y);
        vO = aggregate_win_combined_per_rat(data, sess_types, 'O', ...
            sess_filters{sj}, regions{ri}, resp_only, recruit_thresh, rats_O);
        vY = vY(~isnan(vY)); vO = vO(~isnan(vO));

        if numel(vY) < 2 || numel(vO) < 2
            rows(end+1, :) = {region_titles{ri}, sess_titles{sj}, ...
                numel(vY), numel(vO), NaN, NaN, NaN, NaN, NaN, NaN, NaN};
            continue;
        end
        [~, pval, ~, st] = ttest2(vY, vO);
        rows(end+1, :) = {region_titles{ri}, sess_titles{sj}, ...
            numel(vY), numel(vO), mean(vY), mean(vO), ...
            mean(vY) - mean(vO), NP2_cohens_d(vY, vO, 'unpaired'), ...
            st.tstat, st.df, pval};
    end
end

n_comp = 9;
out = cell2table(rows, 'VariableNames', {'Region','Sess_Filter','N_Young','N_Old', ...
    'Young_Mean','Old_Mean','Mean_Difference','Cohens_d','t_statistic','df','p_value_raw'});
out.p_bonferroni_n9 = min(out.p_value_raw * n_comp, 1);
out.Significant_p05 = out.p_value_raw < 0.05;
end
