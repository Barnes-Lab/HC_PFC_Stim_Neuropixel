function run_full_neuron_analysis_pyr(data, sess_types, ageCols, recruit_thresh, outputDir)
%PYR-only neuron-weighted analysis all
%
%
% INPUT
%   data- per-neuron data (must have .win_combined)
%   sess_types - {'iHC','vHC'}
%   ageCols 
%   recruit_thresh - Hz
%   outputDir 
%
% SS 2026

statsDir = fullfile(outputDir, 'PYR_Stats');
plotsDir = fullfile(outputDir, 'Plots');
if ~exist(statsDir, 'dir'); mkdir(statsDir); end
if ~exist(plotsDir, 'dir'); mkdir(plotsDir); end

T_pyr = NP2_build_neuron_rat_table(data, 'PYR', sess_types, recruit_thresh);
writetable(T_pyr, fullfile(statsDir, 'PYR_rat_level_table.csv'));

NP2_run_neuron_stats(data, sess_types, recruit_thresh, statsDir);

metrics       = {'delta_FR', 'recruitment', 'reliability_per_trial'};
ylabels       = {'\Delta FR (Hz)', '% Recruited', 'Reliability (% trials)'};
within_comps  = {'mono_poly', 'PL_IL', 'iHC_vHC'};
between_grans = {'coarse', 'by_session'};

for m = 1:numel(metrics)
    for c = 1:numel(within_comps)
        plot_within_rat(data, sess_types, within_comps{c}, ...
            metrics{m}, ageCols, plotsDir, ylabels{m}, recruit_thresh);
    end
    for g = 1:numel(between_grans)
        plot_between_age_bars(data, sess_types, between_grans{g}, ...
            metrics{m}, ageCols, plotsDir, ylabels{m}, recruit_thresh);
    end
end

overallStatsDir = fullfile(statsDir, 'overall_FR');
if ~exist(overallStatsDir, 'dir'); mkdir(overallStatsDir); end
modes = {'all', 'responders'};
for mi = 1:2
    plot_overall_FR_per_rat(data, sess_types, modes{mi}, ageCols, plotsDir, recruit_thresh);
    T_stats = NP2_overall_FR_posthoc(data, sess_types, modes{mi}, recruit_thresh);
    writetable(T_stats, fullfile(overallStatsDir, ...
        sprintf('posthoc_overall_FR_%s.csv', modes{mi})));
end

fprintf('  PYR stats + plots + overall-FR: %s\n', outputDir);
end
