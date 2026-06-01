%% ========================================================================
% RUN_NEURON_ANALYSIS_PIPELINE
% ========================================================================
% Spike loading, classification, and response analysis on I80 stims. Stats
% run on PYR only; comparison plots on INT and PYR. Neuron-weighted only.
% 
% SS 2026
% ========================================================================
clear; close all; clc;

%% INPUT PARAMETERS

T_PRM.RATS={'10940','10947','10986','10987','10994',...
    '11039','11042','11043','11040','11044'};
T_PRM.AGES={'O','O','Y','Y','O','O','Y','Y','O','Y'};

T_PRM.BASE_DIR   = 'C:\DATA\NP2';
T_PRM.SR_DIR     = fullfile(T_PRM.BASE_DIR, 'Processed', 'LF_Stimulus_Response_Analysis');
T_PRM.I80_PATH   = fullfile(T_PRM.SR_DIR, 'I80_Analysis_Mono_L5', 'I80_DATA.mat');
T_PRM.STIM_CSV   = fullfile(T_PRM.BASE_DIR, 'ALL_STIM_TIMES_KEEP_CHECK.csv');
T_PRM.MANUAL_CSV = fullfile(T_PRM.SR_DIR, 'manual_window_specifications.csv');
T_PRM.FEPSP_CSV  = fullfile(T_PRM.SR_DIR, 'peak_window_times_post_i80.csv');

T_PRM.WINDOW_MODE    = 'manual';      % 'default' | 'manual' | 'fepsp_peak'
T_PRM.FEPSP_WIN_TYPE = 'tight';       % 'tight' | 'window_all'
T_PRM.WIN_DEFAULTS   = struct('mono', [5 35], 'poly', [35 100]);

if strcmp(T_PRM.WINDOW_MODE, 'fepsp_peak')
    tag = sprintf('fepsp_%s', T_PRM.FEPSP_WIN_TYPE);
else
    tag = T_PRM.WINDOW_MODE;
end
T_PRM.OUT_DIR = fullfile(T_PRM.BASE_DIR, 'Processed', sprintf('Comparisons_Neuron_I80_%s_manual', tag));

T_PRM.PSTH_STIM_MODE = 'I80';         % 'I80' | 'all' | 'both'
T_PRM.BSL_WIN        = [-75, -5];
T_PRM.PL_DEPTH_UM    = 3600;
T_PRM.RECRUIT_THRESH = 1;
T_PRM.APPLY_ARTIFACT = false;

T_PRM.SESS_TYPES   = {'iHC','vHC'};
T_PRM.NEURON_TYPES = {'INT','PYR'};
T_PRM.REGIONS      = {'PL','IL'};

ageCols     = struct('Y', [0.20 0.60 0.20], 'O', [0.50 0.00 0.50]);
classColors = struct('INT', [0.85 0.33 0.10], 'PYR', [0 0.45 0.74]);

if ~exist(T_PRM.OUT_DIR, 'dir'); mkdir(T_PRM.OUT_DIR); end

%% LOAD STIM TIMES + I80

ALL_STIM_TIMES = readtable(T_PRM.STIM_CSV, 'VariableNamingRule', 'preserve');
S_i80    = load(T_PRM.I80_PATH, 'I80_DATA');
I80_DATA = S_i80.I80_DATA;

%% BUILD WINDOW TABLE

window_table = build_neuron_window_table(T_PRM.RATS, T_PRM.SESS_TYPES, T_PRM.REGIONS, ...
    T_PRM.WINDOW_MODE, T_PRM.WIN_DEFAULTS, T_PRM.MANUAL_CSV, ...
    T_PRM.FEPSP_CSV, T_PRM.FEPSP_WIN_TYPE);
writetable(window_table, fullfile(T_PRM.OUT_DIR, ...
    sprintf('window_table_%s.csv', T_PRM.WINDOW_MODE)));

%% LOAD AND CLASSIFY NEURONS

all_SP_filtered        = load_and_filter_spikes(T_PRM.BASE_DIR, T_PRM.RATS, T_PRM.AGES);
neuron_classifications = classify_neurons_filtered(all_SP_filtered, T_PRM.BASE_DIR, T_PRM.APPLY_ARTIFACT);

%% COLLECT NEURON RESPONSES

[data, stim_report] = collect_neuron_data_with_I80(all_SP_filtered, ...
    ALL_STIM_TIMES, neuron_classifications, ...
    window_table, T_PRM.BSL_WIN, T_PRM.PL_DEPTH_UM, ...
    T_PRM.SESS_TYPES, T_PRM.NEURON_TYPES, I80_DATA, T_PRM.RECRUIT_THRESH);
writetable(stim_report, fullfile(T_PRM.OUT_DIR, 'stim_selection_report.csv'));

%% COMPARISON PLOTS

plotsDir = fullfile(T_PRM.OUT_DIR, 'Plots');
if ~exist(plotsDir, 'dir'); mkdir(plotsDir); end
create_combined_window_FR_comparison(data, T_PRM.SESS_TYPES, T_PRM.NEURON_TYPES, ageCols, plotsDir);
create_evoked_neuron_window_comparison(data, T_PRM.SESS_TYPES, T_PRM.NEURON_TYPES, ageCols, plotsDir);

%% PSTH

psthBase = fullfile(T_PRM.OUT_DIR, 'PSTH');
if ~exist(psthBase, 'dir'); mkdir(psthBase); end
run_I80 = any(strcmpi(T_PRM.PSTH_STIM_MODE, {'I80','both'}));
run_all = any(strcmpi(T_PRM.PSTH_STIM_MODE, {'all','both'}));

if run_I80
    T_I80   = filter_stim_times_to_I80(ALL_STIM_TIMES, I80_DATA, T_PRM.SESS_TYPES);
    psthDir = fullfile(psthBase, 'I80');
    if ~exist(psthDir, 'dir'); mkdir(psthDir); end
    create_overall_psth_plots(all_SP_filtered, T_I80, neuron_classifications, ...
        T_PRM.SESS_TYPES, T_PRM.NEURON_TYPES, ageCols, classColors, psthDir);
    create_session_psth_plots(all_SP_filtered, T_I80, neuron_classifications, ...
        T_PRM.SESS_TYPES, T_PRM.NEURON_TYPES, ageCols, classColors, psthDir);
    create_depth_psth_from_comprehensive(all_SP_filtered, T_I80, T_PRM.SESS_TYPES, ageCols, psthDir);
end

if run_all
    psthDir = fullfile(psthBase, 'All');
    if ~exist(psthDir, 'dir'); mkdir(psthDir); end
    create_overall_psth_plots(all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
        T_PRM.SESS_TYPES, T_PRM.NEURON_TYPES, ageCols, classColors, psthDir);
    create_session_psth_plots(all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
        T_PRM.SESS_TYPES, T_PRM.NEURON_TYPES, ageCols, classColors, psthDir);
    create_depth_psth_from_comprehensive(all_SP_filtered, ALL_STIM_TIMES, T_PRM.SESS_TYPES, ageCols, psthDir);
end

%% PER-RAT SANITY CHECK

plot_per_rat_pyr_depth_psth(all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
    I80_DATA, T_PRM.SESS_TYPES, T_PRM.PL_DEPTH_UM, T_PRM.RATS, ...
    fullfile(T_PRM.OUT_DIR, 'PerRat_Sanity'));

%% PYR-ONLY STATS + WITHIN-RAT PLOTS

run_full_neuron_analysis_pyr(data, T_PRM.SESS_TYPES, ageCols, ...
    T_PRM.RECRUIT_THRESH, T_PRM.OUT_DIR);

fprintf('\nDone. Output: %s\n', T_PRM.OUT_DIR);
