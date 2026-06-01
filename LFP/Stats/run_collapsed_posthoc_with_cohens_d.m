function posthoc = run_collapsed_posthoc_with_cohens_d(data, metric_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function posthoc = run_collapsed_posthoc_with_cohens_d(data, metric_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pairwise comparisons using 1 value per rat (n = 5 per age group).
%
% Combined-layer rule (applied via apply_combined_layer_rule when data
% has L2/L5 entries):
%   PL : L5 only
%   IL : mean(L2, L5) per rat
% Sessions are then averaged. Each rat contributes a single value per
% (window x region) cell — matching the per-rat aggregation in bar plots
% and in plot_mono_vs_poly_with_means.
%
% Example call:
%   posthoc = run_collapsed_posthoc_with_cohens_d(data, 'amplitude');
%
% INPUT
%   data        struct      Combined data struct (or already pooled)
%   metric_name char        'amplitude' or 'latency'
%
% OUTPUT
%   posthoc     table       Per-comparison stats with Bonferroni (n=9).
%     Between-age (unpaired): df = 8 (5 young + 5 old)
%     Within-age paired:      df = 4
%     All-rats paired:        df = 9
%
% SS 2026

% Apply the combined-layer rule if data still has the layers field;
% otherwise assume it is already rat-level (e.g. when called by a wrapper
% that pre-pooled).
is_combined = isfield(data, 'layers') && ...
              any(strcmp(data.layers, 'L2')) && ...
              any(strcmp(data.layers, 'L5'));
if is_combined
    data_rat = apply_combined_layer_rule(data);
else
    data_rat = data;
end

if strcmp(metric_name, 'amplitude')
    metric_data = data_rat.amplitude;
else
    metric_data = data_rat.latency;
end

% Average sessions (dim 4) -> [n_rats x window x region]
sess_avg = squeeze(nanmean(metric_data, 4));

mono_PL = sess_avg(:, 1, 1);
mono_IL = sess_avg(:, 1, 2);
poly_PL = sess_avg(:, 2, 1);
poly_IL = sess_avg(:, 2, 2);

% All-channels averages (mean of PL and IL per rat)
mono_all = nanmean([mono_PL mono_IL], 2);
poly_all = nanmean([poly_PL poly_IL], 2);

y_idx = data_rat.young_idx;
o_idx = data_rat.old_idx;

n_comparisons = 9;

C  = {}; G  = {}; N1 = []; N2 = [];
MD = []; CD = []; TS = []; PR = []; PB = []; SIG = {};

% -------- Between-age (Young vs Old) --------
[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_unpaired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    mono_all(y_idx), mono_all(o_idx), 'Mono-AllChannels', 'Young vs Old', n_comparisons);

[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_unpaired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    poly_all(y_idx), poly_all(o_idx), 'Poly-AllChannels', 'Young vs Old', n_comparisons);

[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_unpaired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    poly_PL(y_idx), poly_PL(o_idx), 'Poly-PL', 'Young vs Old', n_comparisons);

[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_unpaired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    poly_IL(y_idx), poly_IL(o_idx), 'Poly-IL', 'Young vs Old', n_comparisons);

% -------- Within-age paired: Mono vs Poly (All Channels) --------
[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_paired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    mono_all, poly_all, 'Mono vs Poly (AllChannels)', 'All', n_comparisons);

[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_paired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    mono_all(y_idx), poly_all(y_idx), 'Mono vs Poly (AllChannels)', 'Young', n_comparisons);

[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_paired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    mono_all(o_idx), poly_all(o_idx), 'Mono vs Poly (AllChannels)', 'Old', n_comparisons);

% -------- Within-age paired: PL vs IL (Poly) --------
[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_paired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    poly_PL, poly_IL, 'PL vs IL (Poly)', 'All', n_comparisons);

[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_paired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    poly_PL(y_idx), poly_IL(y_idx), 'PL vs IL (Poly)', 'Young', n_comparisons);

[C,G,N1,N2,MD,CD,TS,PR,PB,SIG] = add_paired(C,G,N1,N2,MD,CD,TS,PR,PB,SIG, ...
    poly_PL(o_idx), poly_IL(o_idx), 'PL vs IL (Poly)', 'Old', n_comparisons);

posthoc = table(C', G', N1', N2', MD', CD', TS', PR', PB', SIG', ...
    'VariableNames', {'Comparison','Group','N1','N2','Mean_Difference', ...
    'Cohens_d','t_statistic','p_value_raw','p_bonferroni','Significant_Bonf_p05'});

end


%% ========================================================================
function [Comp,Grp,N1,N2,MD,CD,TS,PR,PB,SIG] = add_unpaired( ...
    Comp,Grp,N1,N2,MD,CD,TS,PR,PB,SIG, d1, d2, comp_label, grp_label, n_comp)
v1 = d1(~isnan(d1));
v2 = d2(~isnan(d2));
if length(v1) >= 3 && length(v2) >= 3
    [~, p, ~, stats] = ttest2(v1, v2);
    d = calculate_cohens_d_local(v1, v2, false);
    Comp{end+1} = comp_label;  Grp{end+1}  = grp_label;
    N1(end+1) = length(v1);    N2(end+1) = length(v2);
    MD(end+1) = mean(v1) - mean(v2);  CD(end+1) = d;
    TS(end+1) = stats.tstat;   PR(end+1) = p;
    PB(end+1) = min(p * n_comp, 1);  SIG{end+1} = PB(end) < 0.05;
end
end

%% ========================================================================
function [Comp,Grp,N1,N2,MD,CD,TS,PR,PB,SIG] = add_paired( ...
    Comp,Grp,N1,N2,MD,CD,TS,PR,PB,SIG, d1, d2, comp_label, grp_label, n_comp)
valid = ~isnan(d1) & ~isnan(d2);
v1 = d1(valid);  v2 = d2(valid);
if sum(valid) >= 3
    [~, p, ~, stats] = ttest(v1, v2);
    d = calculate_cohens_d_local(v1, v2, true);
    Comp{end+1} = comp_label;  Grp{end+1}  = grp_label;
    N1(end+1) = sum(valid);    N2(end+1) = sum(valid);
    MD(end+1) = mean(v1 - v2); CD(end+1) = d;
    TS(end+1) = stats.tstat;   PR(end+1) = p;
    PB(end+1) = min(p * n_comp, 1);  SIG{end+1} = PB(end) < 0.05;
end
end

%% ========================================================================
function d = calculate_cohens_d_local(d1, d2, paired)
if paired
    diff = d1(:) - d2(:);
    d = mean(diff) / std(diff);
else
    n1 = length(d1); n2 = length(d2);
    pooled_sd = sqrt(((n1-1)*var(d1) + (n2-1)*var(d2)) / (n1+n2-2));
    d = (mean(d1) - mean(d2)) / pooled_sd;
end
end
