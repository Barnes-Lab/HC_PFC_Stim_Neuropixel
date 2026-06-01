function create_session_psth_plots(all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
    sess_types, neuron_types, ageCols, classColors, figDir)
% CREATE_SESSION_PSTH_PLOTS  For each session (iHC, vHC), two figures:
%   psth_session_<sess>_pooled.{png,fig}      pooled neurons, baseline-subt (Hz)
%   psth_session_<sess>_perRat_dFR.{png,fig}  per-neuron dFR, per-rat avg (Hz)
% SS 2026

time_bins   = -0.500:0.010:0.120;
bin_centers = (time_bins(1:end-1) + time_bins(2:end)) / 2 * 1000;
bsl_idx     = bin_centers >= -400 & bin_centers <= -10;
PL_depth    = 3600;
ages        = {'Y', 'O'};
regions     = {'PL', 'IL'};

for s = 1:length(sess_types)
    sess = sess_types{s};

    data = struct();
    for ti = 1:length(neuron_types)
        nt = neuron_types{ti};
        for ai = 1:2
            a = ages{ai};
            for ri = 1:2
                r = regions{ri};
                data.(nt).(a).(r) = collect_psth_per_cell(all_SP_filtered, ALL_STIM_TIMES, ...
                    neuron_classifications, nt, a, r, {sess}, time_bins, bsl_idx, PL_depth);
            end
        end
    end

    plot_psth_overlay(data, bin_centers, [-10, 120], neuron_types, ageCols, 'pooled', ...
        'Firing Rate Change (Hz)', sprintf('%s Session PSTH (Pooled, biased by N)', sess), ...
        fullfile(figDir, sprintf('psth_session_%s_pooled', sess)));

    plot_psth_overlay(data, bin_centers, [-10, 120], neuron_types, ageCols, 'perRat_bsl', ...
        '\Delta Firing Rate (Hz)', sprintf('%s Session PSTH (Per-Neuron \\Delta FR, Per-Rat Avg)', sess), ...
        fullfile(figDir, sprintf('psth_session_%s_perRat_dFR', sess)));
end
end
