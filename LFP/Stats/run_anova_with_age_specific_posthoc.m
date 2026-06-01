function run_anova_with_age_specific_posthoc(data, T_PRM)
% run_anova_with_age_specific_posthoc(data, T_PRM)
%
% ANOVA + Bonferroni post-hocs 
%
% FOR MONO
%   PL : L5 only
%   IL : mean of L2 and L5 per rat

%   - 3-way ANOVA (Age x Window x Region)            
%   - within-age pairwise     
%   - between-age pairwise    
%   - iHC vs vHC paired (within age)                
%   - between-age per sess        
%
% Layer-specific tests (
%   - L2 vs L5 paired (IL conditions only)           
%   -     
%   -  posthoc_*_poly_collapsed
%
% INPUT:
%   data   struct  amplitude/latency [entry,w,r,s], rat_ids, ages,
%                  young_idx, old_idx; for combined also `layers`.
%   T_PRM  struct  .ANALYSIS_DIR (output directory).
% OUTPUT:
%   CSV files written to T_PRM.ANALYSIS_DIR.
% 2026

fprintf('===Running ANOVA with  Post-hoc Tests ===\n');

is_combined = isfield(data, 'layers') && ...
              any(strcmp(data.layers, 'L2')) && ...
              any(strcmp(data.layers, 'L5'));

if is_combined
    data_pool = apply_combined_layer_rule(data);
else
    data_pool = data;
end

%% Amplitude
write_standard_stats(data_pool, 'amplitude', T_PRM);
if is_combined
    writetable(run_layer_posthoc_with_cohens_d(data, 'amplitude'), ...
        fullfile(T_PRM.ANALYSIS_DIR, 'posthoc_amplitude_layer_paired.csv'));
    writetable(run_between_age_posthoc_by_layer_with_cohens_d(data, 'amplitude'), ...
        fullfile(T_PRM.ANALYSIS_DIR, 'posthoc_amplitude_between_age_by_layer.csv'));
    writetable(run_collapsed_posthoc_with_cohens_d(data, 'amplitude'), ...
        fullfile(T_PRM.ANALYSIS_DIR, 'posthoc_amplitude_poly_collapsed.csv'));
end

%% Latency
write_standard_stats(data_pool, 'latency', T_PRM);
if is_combined
    writetable(run_layer_posthoc_with_cohens_d(data, 'latency'), ...
        fullfile(T_PRM.ANALYSIS_DIR, 'posthoc_latency_layer_paired.csv'));
    writetable(run_between_age_posthoc_by_layer_with_cohens_d(data, 'latency'), ...
        fullfile(T_PRM.ANALYSIS_DIR, 'posthoc_latency_between_age_by_layer.csv'));
    writetable(run_collapsed_posthoc_with_cohens_d(data, 'latency'), ...
        fullfile(T_PRM.ANALYSIS_DIR, 'posthoc_latency_poly_collapsed.csv'));
end

fprintf('Files saved to: %s\n', T_PRM.ANALYSIS_DIR);

end

%% ADDED ALL STATS FUNCTIONS HERE TO BE SPECIFIC AND SO OUTPUT HAS CORRECT LABELS FOR EACH CONDIITON. 
%% WRITE STATS

function write_standard_stats(data_pool, metric, T_PRM)
% Rat-level stat suite (no Layer factor) for one metric.

writetable(run_mixed_anova(data_pool, metric), ...
    fullfile(T_PRM.ANALYSIS_DIR, sprintf('anova_%s_main_effects.csv', metric)));
writetable(run_within_age_posthoc_with_cohens_d(data_pool, metric), ...
    fullfile(T_PRM.ANALYSIS_DIR, sprintf('posthoc_%s_within_age_bonferroni.csv', metric)));
writetable(run_between_age_posthoc_with_cohens_d(data_pool, metric), ...
    fullfile(T_PRM.ANALYSIS_DIR, sprintf('posthoc_%s_between_age_bonferroni.csv', metric)));
writetable(run_session_posthoc_with_cohens_d(data_pool, metric), ...
    fullfile(T_PRM.ANALYSIS_DIR, sprintf('posthoc_%s_session_paired.csv', metric)));
writetable(run_between_age_posthoc_by_session_with_cohens_d(data_pool, metric), ...
    fullfile(T_PRM.ANALYSIS_DIR, sprintf('posthoc_%s_between_age_by_session.csv', metric)));

end


%%  ANOVA ( AGE (between) and window , region within no layer

function anova_results = run_mixed_anova(data, metric_name)
% 3-way mixed ANOVA: Age (between) x Window (within) x Region (within)
% Sessions averaged. One row per rat (or per L2/L5 entry if not pooled).

if strcmp(metric_name, 'amplitude')
    metric_data = data.amplitude;
else
    metric_data = data.latency;
end

avg_data = squeeze(nanmean(metric_data, 4));

Subject = []; Age = {}; Window = []; Region = []; Value = [];
for iEntry = 1:size(avg_data, 1)
    age_val = data.ages{iEntry};
    for w = 1:2
        for r = 1:2
            val = avg_data(iEntry, w, r);
            if ~isnan(val)
                Subject(end+1) = iEntry; %#ok<AGROW>
                Age{end+1}     = age_val;
                Window(end+1)  = w;
                Region(end+1)  = r;
                Value(end+1)   = val;
            end
        end
    end
end

tbl = table(Subject', Age', Window', Region', Value', ...
    'VariableNames', {'Subject', 'Age', 'Window', 'Region', 'Value'});

lme = fitlme(tbl, 'Value ~ Age * Window * Region + (1|Subject)');
anova_results = format_anova_output(anova(lme));

end


%%  Cohen's d 

function d = calculate_cohens_d(data1, data2, paired)
valid1 = ~isnan(data1);
valid2 = ~isnan(data2);

if paired
    valid = valid1 & valid2;
    if sum(valid) < 2, d = NaN; return; end
    diff = data1(valid) - data2(valid);
    d = mean(diff) / std(diff);
else
    d1 = data1(valid1);
    d2 = data2(valid2);
    if length(d1) < 2 || length(d2) < 2, d = NaN; return; end
    n1 = length(d1); n2 = length(d2);
    pooled_sd = sqrt(((n1-1)*var(d1) + (n2-1)*var(d2)) / (n1+n2-2));
    d = (mean(d1) - mean(d2)) / pooled_sd;
end
end



%%  Within-age post-hoc


function posthoc = run_within_age_posthoc_with_cohens_d(data, metric_name)
% Within-age comparisons (Mono vs Poly main effect + 6 pairwise), Bonferroni n=7.

if strcmp(metric_name, 'amplitude')
    metric_data = data.amplitude;
else
    metric_data = data.latency;
end
avg_data = squeeze(nanmean(metric_data, 4));

Age_Group = {}; Comparison = {}; N_Pairs = []; Mean_Diff = [];
Cohens_d = []; T_Stat = []; P_Raw = []; P_Bonf = []; Significant = {};

n_comparisons = 7;

for age_idx = 1:2
    if age_idx == 1
        age_label = 'Young'; entries = data.young_idx;
    else
        age_label = 'Old';   entries = data.old_idx;
    end

    mono_PL = avg_data(entries, 1, 1);
    mono_IL = avg_data(entries, 1, 2);
    poly_PL = avg_data(entries, 2, 1);
    poly_IL = avg_data(entries, 2, 2);

    mono_all = nanmean([mono_PL mono_IL], 2);
    poly_all = nanmean([poly_PL poly_IL], 2);

    [Age_Group, Comparison, N_Pairs, Mean_Diff, Cohens_d, T_Stat, P_Raw, P_Bonf, Significant] = ...
        append_paired(Age_Group, Comparison, N_Pairs, Mean_Diff, Cohens_d, T_Stat, P_Raw, P_Bonf, Significant, ...
        age_label, 'Mono vs Poly (main)', mono_all, poly_all, n_comparisons);

    pairs = {
        mono_PL, mono_IL, 'Mono: PL vs IL';
        mono_PL, poly_PL, 'PL: Mono vs Poly';
        mono_PL, poly_IL, 'Mono-PL vs Poly-IL';
        mono_IL, poly_PL, 'Mono-IL vs Poly-PL';
        mono_IL, poly_IL, 'IL: Mono vs Poly';
        poly_PL, poly_IL, 'Poly: PL vs IL'
    };
    for i = 1:size(pairs, 1)
        [Age_Group, Comparison, N_Pairs, Mean_Diff, Cohens_d, T_Stat, P_Raw, P_Bonf, Significant] = ...
            append_paired(Age_Group, Comparison, N_Pairs, Mean_Diff, Cohens_d, T_Stat, P_Raw, P_Bonf, Significant, ...
            age_label, pairs{i,3}, pairs{i,1}, pairs{i,2}, n_comparisons);
    end
end

posthoc = table(Age_Group', Comparison', N_Pairs', Mean_Diff', Cohens_d', T_Stat', ...
    P_Raw', P_Bonf', Significant', ...
    'VariableNames', {'Age_Group', 'Comparison', 'N_pairs', 'Mean_Difference', ...
    'Cohens_d', 't_statistic', 'p_value_raw', 'p_bonferroni', 'Significant_Bonf_p05'});
end


%%  Between-age post-hoc 

function between_age = run_between_age_posthoc_with_cohens_d(data, metric_name)
fprintf('  Computing between-age post-hoc comparisons with Cohen''s d...\n');

if strcmp(metric_name, 'amplitude')
    metric_data = data.amplitude;
else
    metric_data = data.latency;
end
avg_data = squeeze(nanmean(metric_data, 4));

conditions = {[1 1],'Mono-PL'; [1 2],'Mono-IL'; [2 1],'Poly-PL'; [2 2],'Poly-IL'};
n_comparisons = size(conditions, 1);

Comparison = {}; N_Young = []; N_Old = []; Young_Mean = []; Old_Mean = [];
Mean_Diff = []; Cohens_d = []; T_Stat = []; P_Raw = []; P_Bonf = []; Significant = {};

for i = 1:size(conditions, 1)
    w = conditions{i,1}(1); r = conditions{i,1}(2);
    y = avg_data(data.young_idx, w, r); o = avg_data(data.old_idx, w, r);
    yv = y(~isnan(y)); ov = o(~isnan(o));
    if length(yv) >= 3 && length(ov) >= 3
        [~, p, ~, stats] = ttest2(yv, ov);
        Comparison{end+1} = conditions{i,2};
        N_Young(end+1) = numel(yv); N_Old(end+1) = numel(ov);
        Young_Mean(end+1) = mean(yv); Old_Mean(end+1) = mean(ov);
        Mean_Diff(end+1) = mean(yv) - mean(ov);
        Cohens_d(end+1) = calculate_cohens_d(yv, ov, false);
        T_Stat(end+1) = stats.tstat; P_Raw(end+1) = p;
        P_Bonf(end+1) = min(p * n_comparisons, 1);
        Significant{end+1} = P_Bonf(end) < 0.05;
    end
end

between_age = table(Comparison', N_Young', N_Old', Young_Mean', Old_Mean', ...
    Mean_Diff', Cohens_d', T_Stat', P_Raw', P_Bonf', Significant', ...
    'VariableNames', {'Condition', 'N_Young', 'N_Old', 'Young_Mean', 'Old_Mean', ...
    'Mean_Difference', 'Cohens_d', 't_statistic', 'p_value_raw', 'p_bonferroni', 'Significant_Bonf_p05'});
end


%%   Session post-hoc (iHC vs vHC ) - layer-collapsed 


function session_posthoc = run_session_posthoc_with_cohens_d(data, metric_name)
if strcmp(metric_name, 'amplitude')
    metric_data = data.amplitude;
else
    metric_data = data.latency;
end

conditions = {[1 1],'Mono-PL'; [1 2],'Mono-IL'; [2 1],'Poly-PL'; [2 2],'Poly-IL'};
n_comparisons = size(conditions, 1);
age_groups = {'All', 'Young', 'Old'};

Age_Group = {}; Comparison = {}; N_Pairs = []; iHC_Mean = []; vHC_Mean = [];
Mean_Diff = []; Cohens_d = []; T_Stat = []; P_Raw = []; P_Bonf = []; Significant = {};

for ag = 1:3
    if ag == 1, sel = 1:size(metric_data, 1);
    elseif ag == 2, sel = data.young_idx;
    else, sel = data.old_idx;
    end
    for i = 1:size(conditions, 1)
        w = conditions{i,1}(1); r = conditions{i,1}(2);
        a = metric_data(sel, w, r, 1); b = metric_data(sel, w, r, 2);
        valid = ~isnan(a) & ~isnan(b); av = a(valid); bv = b(valid);
        if sum(valid) >= 3
            [~, p, ~, stats] = ttest(av, bv);
            Age_Group{end+1} = age_groups{ag};
            Comparison{end+1} = conditions{i,2};
            N_Pairs(end+1) = sum(valid);
            iHC_Mean(end+1) = mean(av); vHC_Mean(end+1) = mean(bv);
            Mean_Diff(end+1) = mean(av - bv);
            Cohens_d(end+1) = calculate_cohens_d(av, bv, true);
            T_Stat(end+1) = stats.tstat; P_Raw(end+1) = p;
            P_Bonf(end+1) = min(p * n_comparisons, 1);
            Significant{end+1} = P_Bonf(end) < 0.05;
        end
    end
end

session_posthoc = table(Age_Group', Comparison', N_Pairs', iHC_Mean', vHC_Mean', ...
    Mean_Diff', Cohens_d', T_Stat', P_Raw', P_Bonf', Significant', ...
    'VariableNames', {'Age_Group', 'Condition', 'N_pairs', 'iHC_Mean', 'vHC_Mean', ...
    'Mean_Difference', 'Cohens_d', 't_statistic', 'p_value_raw', 'p_bonferroni', 'Significant_Bonf_p05'});
end


%%  Layer post-hoc (L2 vs L5 ) - IL only

function layer_posthoc = run_layer_posthoc_with_cohens_d(data, metric_name)
% Per combined-layer rule, PL is L5-only (no L2 vs L5 contrast).
% Restricted to IL conditions: Mono-IL, Poly-IL. Split by age: All/Young/Old.

if strcmp(metric_name, 'amplitude')
    metric_data = data.amplitude;
else
    metric_data = data.latency;
end
avg_data = squeeze(nanmean(metric_data, 4));

unique_rats = unique(data.rat_ids);
n_rats = numel(unique_rats);
rat_age = cell(n_rats, 1);
for iRat = 1:n_rats
    any_idx = find(strcmp(data.rat_ids, unique_rats{iRat}), 1);
    rat_age{iRat} = data.ages{any_idx};
end

conditions = {[1 2],'Mono-IL'; [2 2],'Poly-IL'};
n_comparisons = size(conditions, 1);
age_groups = {'All', 'Young', 'Old'};

Age_Group = {}; Comparison = {}; N_Pairs = []; L2_Mean = []; L5_Mean = [];
Mean_Diff = []; Cohens_d = []; T_Stat = []; P_Raw = []; P_Bonf = []; Significant = {};

for ag = 1:3
    if ag == 1, rat_sel = true(n_rats, 1);
    elseif ag == 2, rat_sel = strcmp(rat_age, 'Y');
    else, rat_sel = strcmp(rat_age, 'O');
    end

    for i = 1:size(conditions, 1)
        w = conditions{i,1}(1); r = conditions{i,1}(2);

        L2_vals = nan(n_rats, 1); L5_vals = nan(n_rats, 1);
        for iRat = 1:n_rats
            rat = unique_rats{iRat};
            L2_idx = find(strcmp(data.rat_ids, rat) & strcmp(data.layers, 'L2'), 1);
            L5_idx = find(strcmp(data.rat_ids, rat) & strcmp(data.layers, 'L5'), 1);
            if ~isempty(L2_idx) && ~isempty(L5_idx)
                L2_vals(iRat) = avg_data(L2_idx, w, r);
                L5_vals(iRat) = avg_data(L5_idx, w, r);
            end
        end
        L2_vals = L2_vals(rat_sel); L5_vals = L5_vals(rat_sel);

        valid = ~isnan(L2_vals) & ~isnan(L5_vals);
        a = L2_vals(valid); b = L5_vals(valid);
        if sum(valid) >= 3
            [~, p, ~, stats] = ttest(a, b);
            Age_Group{end+1} = age_groups{ag};
            Comparison{end+1} = conditions{i,2};
            N_Pairs(end+1) = sum(valid);
            L2_Mean(end+1) = mean(a); L5_Mean(end+1) = mean(b);
            Mean_Diff(end+1) = mean(a - b);
            Cohens_d(end+1) = calculate_cohens_d(a, b, true);
            T_Stat(end+1) = stats.tstat; P_Raw(end+1) = p;
            P_Bonf(end+1) = min(p * n_comparisons, 1);
            Significant{end+1} = P_Bonf(end) < 0.05;
        end
    end
end

layer_posthoc = table(Age_Group', Comparison', N_Pairs', L2_Mean', L5_Mean', ...
    Mean_Diff', Cohens_d', T_Stat', P_Raw', P_Bonf', Significant', ...
    'VariableNames', {'Age_Group', 'Condition', 'N_pairs', 'L2_Mean', 'L5_Mean', ...
    'Mean_Difference', 'Cohens_d', 't_statistic', 'p_value_raw', 'p_bonferroni', 'Significant_Bonf_p05'});
end


%%  Between-age by Layer - L2- IL only, L5-all conditions

function between_age_layer = run_between_age_posthoc_by_layer_with_cohens_d(data, metric_name)
% Young vs Old stratified by Layer. Per combined-layer rule:
%   L2: IL conditions only (Mono-IL, Poly-IL)
%   L5: all conditions (Mono-PL, Mono-IL, Poly-PL, Poly-IL)
% Bonferroni applied per-layer (within its own family).

if strcmp(metric_name, 'amplitude')
    metric_data = data.amplitude;
else
    metric_data = data.latency;
end
avg_data = squeeze(nanmean(metric_data, 4));

all_conditions = {[1 1],'Mono-PL'; [1 2],'Mono-IL'; [2 1],'Poly-PL'; [2 2],'Poly-IL'};
IL_conditions  = {[1 2],'Mono-IL'; [2 2],'Poly-IL'};

Layer = {}; Comparison = {}; N_Young = []; N_Old = []; Young_Mean = []; Old_Mean = [];
Mean_Diff = []; Cohens_d = []; T_Stat = []; P_Raw = []; P_Bonf = []; Significant = {};

layers = {'L2', 'L5'};
for L = 1:2
    layer_label = layers{L};
    if strcmp(layer_label, 'L2')
        conditions = IL_conditions;
    else
        conditions = all_conditions;
    end
    n_comparisons = size(conditions, 1);

    layer_mask = strcmp(data.layers, layer_label);
    young_sel = intersect(data.young_idx, find(layer_mask));
    old_sel   = intersect(data.old_idx,   find(layer_mask));

    for i = 1:size(conditions, 1)
        w = conditions{i,1}(1); r = conditions{i,1}(2);
        y = avg_data(young_sel, w, r); o = avg_data(old_sel, w, r);
        yv = y(~isnan(y)); ov = o(~isnan(o));
        if length(yv) >= 3 && length(ov) >= 3
            [~, p, ~, stats] = ttest2(yv, ov);
            Layer{end+1} = layer_label;
            Comparison{end+1} = conditions{i,2};
            N_Young(end+1) = numel(yv); N_Old(end+1) = numel(ov);
            Young_Mean(end+1) = mean(yv); Old_Mean(end+1) = mean(ov);
            Mean_Diff(end+1) = mean(yv) - mean(ov);
            Cohens_d(end+1) = calculate_cohens_d(yv, ov, false);
            T_Stat(end+1) = stats.tstat; P_Raw(end+1) = p;
            P_Bonf(end+1) = min(p * n_comparisons, 1);
            Significant{end+1} = P_Bonf(end) < 0.05;
        end
    end
end

between_age_layer = table(Layer', Comparison', N_Young', N_Old', Young_Mean', Old_Mean', ...
    Mean_Diff', Cohens_d', T_Stat', P_Raw', P_Bonf', Significant', ...
    'VariableNames', {'Layer', 'Condition', 'N_Young', 'N_Old', 'Young_Mean', 'Old_Mean', ...
    'Mean_Difference', 'Cohens_d', 't_statistic', 'p_value_raw', 'p_bonferroni', 'Significant_Bonf_p05'});
end


%%  Between-age posthoc- per sess


function between_age_session = run_between_age_posthoc_by_session_with_cohens_d(data, metric_name)
% Receives rat-level pooled data (combined rule already applied upstream).

is_combined = isfield(data, 'layers') && ...
              any(strcmp(data.layers, 'L2')) && ...
              any(strcmp(data.layers, 'L5'));
if is_combined
    data = apply_combined_layer_rule(data);
end

if strcmp(metric_name, 'amplitude')
    metric_data = data.amplitude;
else
    metric_data = data.latency;
end

conditions = {[1 1],'Mono-PL'; [1 2],'Mono-IL'; [2 1],'Poly-PL'; [2 2],'Poly-IL'};
n_comparisons = size(conditions, 1);
sessions = {'iHC', 'vHC'};

Session = {}; Comparison = {}; N_Young = []; N_Old = []; Young_Mean = []; Old_Mean = [];
Mean_Diff = []; Cohens_d = []; T_Stat = []; P_Raw = []; P_Bonf = []; Significant = {};

for s = 1:2
    for i = 1:size(conditions, 1)
        w = conditions{i,1}(1); r = conditions{i,1}(2);
        y = metric_data(data.young_idx, w, r, s);
        o = metric_data(data.old_idx,   w, r, s);
        yv = y(~isnan(y)); ov = o(~isnan(o));
        if length(yv) >= 3 && length(ov) >= 3
            [~, p, ~, stats] = ttest2(yv, ov);
            Session{end+1} = sessions{s};
            Comparison{end+1} = conditions{i,2};
            N_Young(end+1) = numel(yv); N_Old(end+1) = numel(ov);
            Young_Mean(end+1) = mean(yv); Old_Mean(end+1) = mean(ov);
            Mean_Diff(end+1) = mean(yv) - mean(ov);
            Cohens_d(end+1) = calculate_cohens_d(yv, ov, false);
            T_Stat(end+1) = stats.tstat; P_Raw(end+1) = p;
            P_Bonf(end+1) = min(p * n_comparisons, 1);
            Significant{end+1} = P_Bonf(end) < 0.05;
        end
    end
end

between_age_session = table(Session', Comparison', N_Young', N_Old', Young_Mean', Old_Mean', ...
    Mean_Diff', Cohens_d', T_Stat', P_Raw', P_Bonf', Significant', ...
    'VariableNames', {'Session', 'Condition', 'N_Young', 'N_Old', 'Young_Mean', 'Old_Mean', ...
    'Mean_Difference', 'Cohens_d', 't_statistic', 'p_value_raw', 'p_bonferroni', 'Significant_Bonf_p05'});
end


%%  OTHERS
function tbl = format_anova_output(anova_output)
if isa(anova_output, 'dataset')
    Effect = cellstr(anova_output.Term);
    F_Stat = double(anova_output.FStat);
    DF1 = double(anova_output.DF1); DF2 = double(anova_output.DF2);
    P_Value = double(anova_output.pValue);
elseif istable(anova_output)
    Effect = anova_output.Properties.RowNames;
    F_Stat = anova_output.FStat;
    DF1 = anova_output.DF1; DF2 = anova_output.DF2;
    P_Value = anova_output.pValue;
else
    error('Unknown ANOVA output format');
end
Significant = P_Value < 0.05;
tbl = table(Effect, F_Stat, DF1, DF2, P_Value, Significant, ...
    'VariableNames', {'Effect', 'F_statistic', 'df1', 'df2', 'p_value', 'Significant_p05'});
end

%% SAVE ALL
function [Age_Group, Comparison, N_Pairs, Mean_Diff, Cohens_d, T_Stat, P_Raw, P_Bonf, Significant] = ...
    append_paired(Age_Group, Comparison, N_Pairs, Mean_Diff, Cohens_d, T_Stat, P_Raw, P_Bonf, Significant, ...
    age_label, comp_name, d1, d2, n_comparisons)
% Appends a paired-test row (d2 - d1) to the within-age post-hoc accumulators.

valid = ~isnan(d1) & ~isnan(d2);
if sum(valid) < 3, return; end
[~, p, ~, stats] = ttest(d2(valid), d1(valid));
Age_Group{end+1}   = age_label;
Comparison{end+1}  = comp_name;
N_Pairs(end+1)     = sum(valid);
Mean_Diff(end+1)   = mean(d2(valid) - d1(valid));
Cohens_d(end+1)    = calculate_cohens_d(d2, d1, true);
T_Stat(end+1)      = stats.tstat;
P_Raw(end+1)       = p;
P_Bonf(end+1)      = min(p * n_comparisons, 1);
Significant{end+1} = P_Bonf(end) < 0.05;
end
