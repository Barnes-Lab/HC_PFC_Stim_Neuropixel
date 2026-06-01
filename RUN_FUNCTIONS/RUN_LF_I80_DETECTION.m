%% ========================================================================
% RUN_LF_I80_DETECTION
% ========================================================================
% I80 detection on IL mono amplitudes. Per-shank fits, hierarchy selection,
% flat I80 table + per-rat verification PNGs + age comparison plots.
% SS 2026
% ========================================================================
clear; close all; clc;

%% INPUT PARAMETERS

T_PRM.RATS={'10940','10947','10986','10987','10994',...
    '11039','11042','11043','11040','11044'};
T_PRM.AGES={'O','O','Y','Y','O','O','Y','Y','O','Y'};

T_PRM.STORAGE_DIR    = 'C:\DATA\NP2\Processed\LF_Stimulus_Response_Analysis';
T_PRM.MONO_DATA_PATH = fullfile(T_PRM.STORAGE_DIR, 'MONO_DATA.mat');
T_PRM.STIM_CSV_PATH  = fullfile(T_PRM.STORAGE_DIR, '..', 'ALL_STIM_TIMES_KEEP_CHECK.csv');
T_PRM.OUTPUT_DIR     = fullfile(T_PRM.STORAGE_DIR, 'I80_Analysis_Mono_L5');
T_PRM.VERIFY_DIR     = fullfile(T_PRM.OUTPUT_DIR, 'verify');

T_PRM.PRM_R2_min             = 0.40;
T_PRM.PRM_R2_forced          = 0.10;
T_PRM.PRM_I80_range          = [50, 900];
T_PRM.PRM_n_analysis         = 5;
T_PRM.PRM_flat_thresh_mV     = 0.10;
T_PRM.PRM_flat_force_uA      = 200;
T_PRM.PRM_I80_diff_thresh_uA = 100;
T_PRM.PRM_amp_field          = 'mean_amplitude';

PARAMS.colors.young = [0, 0.6, 0.3];
PARAMS.colors.old   = [0.5, 0, 0.5];

%% PER-RAT OVERRIDES

T_PRM.OVR = struct();
T_PRM.OVR.Rat10940.iHC = struct('shank', 4, 'force_method', 'auto');
T_PRM.OVR.Rat10940.vHC = struct('shank', 4, 'force_method', 'auto');
T_PRM.OVR.Rat10947.iHC = struct('shank', 3, 'force_method', 'auto');
T_PRM.OVR.Rat10947.vHC = struct('shank', 3, 'force_method', 'exp_sat');
T_PRM.OVR.Rat10986.iHC = struct('shank', 3, 'force_method', 'auto');
T_PRM.OVR.Rat10986.vHC = struct('shank', 3, 'force_method', 'auto');
T_PRM.OVR.Rat10987.iHC = struct('shank', 3, 'force_method', 'auto');
T_PRM.OVR.Rat10987.vHC = struct('shank', 3, 'force_method', 'auto');
T_PRM.OVR.Rat10994.iHC = struct('shank', 4, 'force_method', 'auto');
T_PRM.OVR.Rat10994.vHC = struct('shank', 4, 'force_method', 'auto');
T_PRM.OVR.Rat10995.iHC = struct('shank', 3, 'force_method', 'exp_sat');
T_PRM.OVR.Rat10995.vHC = struct('shank', 2, 'force_method', 'auto');
T_PRM.OVR.Rat11039.iHC = struct('shank', 4, 'force_method', 'auto');
T_PRM.OVR.Rat11039.vHC = struct('shank', 4, 'force_method', 'auto');
T_PRM.OVR.Rat11040.iHC = struct('shank', 3, 'force_method', 'auto');
T_PRM.OVR.Rat11040.vHC = struct('shank', 3, 'force_method', 'exp_sat');
T_PRM.OVR.Rat11042.iHC = struct('shank', 4, 'force_method', 'auto');
T_PRM.OVR.Rat11042.vHC = struct('shank', 4, 'force_method', 'auto');
T_PRM.OVR.Rat11043.iHC = struct('shank', 2, 'force_method', 'auto');
T_PRM.OVR.Rat11043.vHC = struct('shank', 2, 'force_flat', true);
T_PRM.OVR.Rat11044.iHC = struct('shank', 2, 'force_method', 'auto');
T_PRM.OVR.Rat11044.vHC = struct('shank', 2, 'force_method', 'auto');

if ~exist(T_PRM.OUTPUT_DIR, 'dir'); mkdir(T_PRM.OUTPUT_DIR); end
if ~exist(T_PRM.VERIFY_DIR, 'dir'); mkdir(T_PRM.VERIFY_DIR); end

%% LOAD

load(T_PRM.MONO_DATA_PATH, 'MONO_DATA');
stim_tbl = readtable(T_PRM.STIM_CSV_PATH, 'VariableNamingRule', 'preserve');

%% FIT PER SHANK + DECISION

I80_DATA = struct();
for iRat = 1:numel(T_PRM.RATS)
    rat_id = T_PRM.RATS{iRat};
    age    = T_PRM.AGES{iRat};
    ratFld = sprintf('Rat%s', rat_id);
    if ~isfield(MONO_DATA, ratFld); continue; end

    L5_shank = NaN; L2_shank = NaN;
    if isfield(MONO_DATA.(ratFld), 'L5_shank'); L5_shank = MONO_DATA.(ratFld).L5_shank; end
    if isfield(MONO_DATA.(ratFld), 'L2_shank'); L2_shank = MONO_DATA.(ratFld).L2_shank; end

    for sess_idx = 1:2
        sess_name = peek_session_name(MONO_DATA.(ratFld), sess_idx);
        if isempty(sess_name); continue; end
        ovr = get_override(T_PRM.OVR, ratFld, sess_name);

        shank_results = repmat(empty_shank_result(), 4, 1);
        for sh = 1:4
            shFld = sprintf('Shank%d', sh);
            if ~isfield(MONO_DATA.(ratFld), shFld); continue; end
            force_method = '';
            if ~isempty(ovr) && isfield(ovr, 'shank') && ovr.shank == sh && ...
                    isfield(ovr, 'force_method')
                force_method = ovr.force_method;
            end
            shank_results(sh) = fit_one_shank_i80( ...
                MONO_DATA.(ratFld).(shFld), sh, sess_idx, stim_tbl, rat_id, ...
                'PRM_R2_min',         T_PRM.PRM_R2_min, ...
                'PRM_R2_forced',      T_PRM.PRM_R2_forced, ...
                'PRM_I80_range',      T_PRM.PRM_I80_range, ...
                'PRM_flat_thresh_mV', T_PRM.PRM_flat_thresh_mV, ...
                'PRM_amp_field',      T_PRM.PRM_amp_field, ...
                'PRM_force_method',   force_method);
        end

        decision = apply_shank_hierarchy(shank_results, L5_shank, ovr, ...
            'PRM_I80_diff_thresh_uA', T_PRM.PRM_I80_diff_thresh_uA, ...
            'PRM_flat_force_uA',      T_PRM.PRM_flat_force_uA);

        as = struct('stim_idx_valid', [], 'stim_idx_sorted', [], 'currents', []);
        if ~isnan(decision.chosen_shank) && ~isnan(decision.I80)
            sr = shank_results(decision.chosen_shank);
            if ~isempty(sr.stim_currents_all) && ~isempty(sr.valid_stims)
                as = select_analysis_stims(decision.I80, sr.stim_currents_all, ...
                    sr.valid_stims, T_PRM.PRM_n_analysis);
            end
        end

        I80_DATA.(ratFld).(sess_name).age                = age;
        I80_DATA.(ratFld).(sess_name).shank_results      = shank_results;
        I80_DATA.(ratFld).(sess_name).chosen_shank       = decision.chosen_shank;
        I80_DATA.(ratFld).(sess_name).chosen_source      = decision.chosen_source;
        I80_DATA.(ratFld).(sess_name).I80                = decision.I80;
        I80_DATA.(ratFld).(sess_name).L5_I80             = decision.L5_I80;
        I80_DATA.(ratFld).(sess_name).Best_R2_Shank      = decision.Best_R2_Shank;
        I80_DATA.(ratFld).(sess_name).L5_vs_Best_diff_uA = decision.L5_vs_Best_diff_uA;
        I80_DATA.(ratFld).(sess_name).analysis_stims     = as;
    end
end

%% SAVE TABLE

I80_TABLE = build_i80_table_flat(I80_DATA, T_PRM.RATS, T_PRM.AGES);
save(fullfile(T_PRM.OUTPUT_DIR, 'I80_DATA.mat'), 'I80_DATA', 'I80_TABLE', 'T_PRM');
writetable(I80_TABLE, fullfile(T_PRM.OUTPUT_DIR, 'I80_TABLE.csv'));
fprintf('Saved I80_DATA.mat + I80_TABLE.csv\n');

%% VERIFICATION FIGURES

for iRat = 1:numel(T_PRM.RATS)
    rat_id = T_PRM.RATS{iRat};
    ratFld = sprintf('Rat%s', rat_id);
    if ~isfield(I80_DATA, ratFld); continue; end
    sess_names = fieldnames(I80_DATA.(ratFld));
    for s = 1:numel(sess_names)
        sn  = sess_names{s};
        d   = I80_DATA.(ratFld).(sn);
        png = fullfile(T_PRM.VERIFY_DIR, sprintf('Rat%s_%s_I80_verify.png', rat_id, sn));
        plot_i80_verification(rat_id, sn, d.shank_results, ...
            d.chosen_shank, d.I80, png, ...
            'PRM_chosen_source', d.chosen_source);
    end
end

%% AGE COMPARISON FIGURES

plot_i80_paper_perrat(I80_DATA, T_PRM, T_PRM.OUTPUT_DIR);
plot_I80_scatter_L5(I80_TABLE, PARAMS, T_PRM.OUTPUT_DIR, T_PRM.RATS);
plot_age_comparison_normalized(I80_TABLE, MONO_DATA, PARAMS, ...
    T_PRM.OUTPUT_DIR, T_PRM.RATS, T_PRM.PRM_amp_field);
