function out = NP2_neuron_anova(data, sess_types, metric, recruit_thresh)
% NP2_NEURON_ANOVA  3-way rm ANOVA Age x Window x Region with random
%
% SS 2026

regions = {'PL', 'IL'};
windows = {'mono', 'poly'};
ages    = {'Y', 'O'};

rows = {};
for ai = 1:2
    a = ages{ai};
    rats_list = unique_rats_for_age(data, sess_types, 'PYR', a);
    for ri = 1:2
        for wi = 1:2
            v = get_per_rat_value(data, sess_types, a, 'all', regions{ri}, windows{wi}, ...
                metric, recruit_thresh, rats_list);
            for r = 1:numel(rats_list)
                if isnan(v(r)); continue; end
                rows(end+1, :) = {rats_list{r}, a, regions{ri}, windows{wi}, v(r)};
            end
        end
    end
end

T_avg = cell2table(rows, 'VariableNames', {'Rat','Age','Region','Window','Y'});
T_avg.Age    = categorical(T_avg.Age,    {'Y','O'});
T_avg.Region = categorical(T_avg.Region, {'PL','IL'});
T_avg.Window = categorical(T_avg.Window, {'mono','poly'});
T_avg.Rat    = categorical(T_avg.Rat);

mdl = fitlme(T_avg, 'Y ~ Age * Window * Region + (1|Rat)', 'DummyVarCoding', 'effects');
A   = anova(mdl);
eff  = string(A.Term);
keep = eff ~= "(Intercept)";
out  = table(cellstr(eff(keep)), A.FStat(keep), A.DF1(keep), A.DF2(keep), ...
    A.pValue(keep), A.pValue(keep) < 0.05, ...
    'VariableNames', {'Effect','F_statistic','df1','df2','p_value','Significant_p05'});
end
