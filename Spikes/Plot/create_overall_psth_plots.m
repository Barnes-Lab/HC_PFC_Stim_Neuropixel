function create_overall_psth_plots(all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
    sess_types, neuron_types, ageCols, classColors, figDir)
% CREATE_OVERALL_PSTH_PLOTS  Two figures, sessions pooled:
%   psth_overall_pooled.{png,fig}        baseline-subtracted (Hz), n=neurons (biased)
%   psth_overall_perRat_dFR.{png,fig}    per-neuron dFR, per-rat avg (Hz, unbiased)
% Baseline z-score variant deprecated (was inflating cross-rat variance).
% SS 2026

time_bins   = -0.500:0.010:0.120;
bin_centers = (time_bins(1:end-1) + time_bins(2:end)) / 2 * 1000;
bsl_idx     = bin_centers >= -400 & bin_centers <= -10;
PL_depth    = 3600;
ages        = {'Y', 'O'};
regions     = {'PL', 'IL'};

data = struct();
for ti = 1:length(neuron_types)
    nt = neuron_types{ti};
    for ai = 1:2
        a = ages{ai};
        for ri = 1:2
            r = regions{ri};
            data.(nt).(a).(r) = collect_psth_per_cell(all_SP_filtered, ALL_STIM_TIMES, ...
                neuron_classifications, nt, a, r, sess_types, time_bins, bsl_idx, PL_depth);
        end
    end
end

plot_psth_overlay(data, bin_centers, [-20, 120], neuron_types, ageCols, 'pooled', ...
    'Firing Rate Change (Hz)', 'Overall PSTH (Pooled Neurons, biased by N)', ...
    fullfile(figDir, 'psth_overall_pooled'));

plot_psth_overlay(data, bin_centers, [-20, 120], neuron_types, ageCols, 'perRat_bsl', ...
    '\Delta Firing Rate (Hz)', 'Overall PSTH (Per-Neuron \Delta FR, Per-Rat Avg)', ...
    fullfile(figDir, 'psth_overall_perRat_dFR'));
end
