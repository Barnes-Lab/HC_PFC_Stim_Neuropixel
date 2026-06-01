%% ========================================================================
% RUN_LF_COMPARISONS_PLOT_STATS
% ========================================================================
% Comparisons of fEPSP properties across 5 I80-selected stims per recording.
% Runs Combined-layer condition; generates figures and statistics.
% SS 2026
% ========================================================================
clear; close all; clc;

%% INPUT PARAMETERS

T_PRM.RATS={'10940','10947','10986','10987','10994',...
    '11039','11042','11043','11040','11044'};
T_PRM.AGES={'O','O','Y','Y','O','O','Y','Y','O','Y'};

T_PRM.n_young = sum(strcmp(T_PRM.AGES, 'Y'));
T_PRM.n_old   = sum(strcmp(T_PRM.AGES, 'O'));

T_PRM.STORAGE_DIR     = 'C:\DATA\NP2\Processed';
T_PRM.SR_ANALYSIS_DIR = fullfile(T_PRM.STORAGE_DIR, 'LF_Stimulus_Response_Analysis');
T_PRM.IO_ANALYSIS_DIR = fullfile(T_PRM.SR_ANALYSIS_DIR, 'I80_Analysis_Mono_L5');
T_PRM.BASE_OUTPUT_DIR = fullfile(T_PRM.SR_ANALYSIS_DIR, '35ms_cutoff_COMP_FIGS_lowest_v2/');

amplitude_var = 'lowest';   % 'lowest' | '3chan' | '50pct' | '25pct'

if ~exist(T_PRM.BASE_OUTPUT_DIR, 'dir'); mkdir(T_PRM.BASE_OUTPUT_DIR); end

COLORS = struct();
COLORS.young = [0, 0.6, 0.3];
COLORS.old   = [0.5, 0, 0.5];
COLORS.PL    = [0.3, 0.7, 0.9];
COLORS.IL    = [0.9, 0.5, 0.2];
COLORS.mono  = [0.2, 0.5, 0.2];
COLORS.poly  = [0.6, 0.3, 0.6];

%% LOAD DATA

sr_file = fullfile(T_PRM.SR_ANALYSIS_DIR, 'SR_DATA_peaks.mat');
if ~exist(sr_file, 'file')
    error('SR_DATA_peaks.mat not found. Run RUN_LF_SLOPE_DETECTION first.');
end
load(sr_file, 'SR_DATA', 'PARAMS');

%% RUN PER LAYER CONDITION

layer_conditions = {'Combined'};
for layer_idx = 1:length(layer_conditions)
    layer_filter = layer_conditions{layer_idx};
    fprintf('  ANALYZING: %s\n', layer_filter);

    T_PRM.ANALYSIS_DIR = fullfile(T_PRM.BASE_OUTPUT_DIR, [layer_filter '_Analysis']);
    T_PRM.FIG_DIR      = fullfile(T_PRM.ANALYSIS_DIR, 'Figures');
    T_PRM.PUB_FIG_DIR  = T_PRM.FIG_DIR;
    if ~exist(T_PRM.ANALYSIS_DIR, 'dir'); mkdir(T_PRM.ANALYSIS_DIR); end
    if ~exist(T_PRM.FIG_DIR,      'dir'); mkdir(T_PRM.FIG_DIR);      end

    [data, stim_report] = extract_negative_data(SR_DATA, T_PRM, amplitude_var, layer_filter);

    run_anova_with_age_specific_posthoc(data, T_PRM);
    create_all_I80_figures(data, COLORS, T_PRM);

    fprintf('\n=== %s Analysis Complete ===\n', layer_filter);
    fprintf('Output: %s\n', T_PRM.ANALYSIS_DIR);
end
