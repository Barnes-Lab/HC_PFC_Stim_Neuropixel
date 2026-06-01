function NP2_run_neuron_stats(data, sess_types, recruit_thresh, outDir)
% NP2_RUN_NEURON_STATS  
%   <outDir>/<metric>/
%       anova_<metric>_main_effects.csv
%       posthoc_<metric>_within_age.csv
%       posthoc_<metric>_between_age.csv
%       posthoc_<metric>_session_paired.csv
%       posthoc_<metric>_between_age_by_session.csv
%
% INPUT
%   data -- per-neuron data from collect_neuron_data_with_I80
%   sess_types
%   recruit_thresh - Hz; threshold for recruitment / reliability
%   outDir -PYR_Stats root
%
% SS 2026

metrics = {'delta_FR', 'recruitment', 'reliability_per_trial'};
if ~exist(outDir, 'dir'); mkdir(outDir); end

for m = 1:numel(metrics)
    metric = metrics{m};
    mDir   = fullfile(outDir, metric);
    if ~exist(mDir, 'dir'); mkdir(mDir); end

    A  = NP2_neuron_anova                         (data, sess_types, metric, recruit_thresh);
    W  = NP2_neuron_within_age_posthoc            (data, sess_types, metric, recruit_thresh);
    B  = NP2_neuron_between_age_posthoc           (data, sess_types, metric, recruit_thresh);
    S  = NP2_neuron_session_paired_posthoc        (data, sess_types, metric, recruit_thresh);
    BS = NP2_neuron_between_age_by_session_posthoc(data, sess_types, metric, recruit_thresh);

    writetable(A,  fullfile(mDir, sprintf('anova_%s_main_effects.csv',                metric)));
    writetable(W,  fullfile(mDir, sprintf('posthoc_%s_within_age.csv',                metric)));
    writetable(B,  fullfile(mDir, sprintf('posthoc_%s_between_age.csv',               metric)));
    writetable(S,  fullfile(mDir, sprintf('posthoc_%s_session_paired.csv',            metric)));
    writetable(BS, fullfile(mDir, sprintf('posthoc_%s_between_age_by_session.csv',    metric)));
end
end
