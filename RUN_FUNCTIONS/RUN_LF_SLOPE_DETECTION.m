%% ========================================================================
% RUN_LF_SLOPE_DETECTION
% ========================================================================
% Negative peak detection in MONO and POLY windows on PL and IL channels
% for L2 and L5 shanks. Runs on the 5 I80-selected analysis stims.
% SS 2026
% ========================================================================
clear; close all; clc;

%% INPUT PARAMETERS


T_PRM.RATS={'10940','10947','10986','10987','10994',...
    '11039','11042','11043','11040','11044'};
T_PRM.AGES={'O','O','Y','Y','O','O','Y','Y','O','Y'};

T_PRM.STORAGE_DIR   = 'C:\DATA\NP2\Processed';
T_PRM.OUTPUT_DIR    = fullfile(T_PRM.STORAGE_DIR, 'LF_Stimulus_Response_Analysis');
T_PRM.I80_DIR       = fullfile(T_PRM.OUTPUT_DIR, 'I80_Analysis_Mono_L5');
T_PRM.FIG_DIR       = fullfile(T_PRM.OUTPUT_DIR, 'SlopeDetection_AfterI80_35ms');
T_PRM.STIM_CSV      = 'ALL_STIM_TIMES_KEEP_CHECK.csv';
T_PRM.WINDOW_CSV    = 'manual_window_specifications.csv';
T_PRM.SHANK_MAP_CSV = fullfile(T_PRM.STORAGE_DIR, 'ShankInfo', 'ShankMap.csv');
T_PRM.CHAN_INCL_DIR = fullfile(T_PRM.OUTPUT_DIR, 'ChannelInclusion');

PARAMS = struct();
PARAMS.IL_mono_window_ms          = [10, 35];
PARAMS.IL_poly_window_ms          = [35, 100];
PARAMS.PL_mono_window_ms          = [10, 45];
PARAMS.PL_poly_window_ms          = [45, 110];
PARAMS.PL_IL_boundary_um          = 3600;
PARAMS.baseline_window_ms         = [-205, -5];
PARAMS.bad_stim_thresh            = 0.35;
PARAMS.prepeak_baseline_ms        = [5, 8];
PARAMS.smooth_window_ms           = 3;
PARAMS.drift_thresh_mV            = 0.3;
PARAMS.amp_frac_thresh            = 0.5;
PARAMS.min_peak_prominence_mV     = 0.05;
PARAMS.per_trial_half_win_mono_ms = 5;
PARAMS.per_trial_half_win_poly_ms = 20;
PARAMS.colors.PL                  = [0.3, 0.7, 0.9];
PARAMS.colors.IL                  = [0.88, 0.35, 0.11];

if ~exist(T_PRM.OUTPUT_DIR, 'dir'); mkdir(T_PRM.OUTPUT_DIR); end
if ~exist(T_PRM.FIG_DIR,    'dir'); mkdir(T_PRM.FIG_DIR);    end

%% LOAD AUXILIARY TABLES

allStimTbl = readtable(fullfile(T_PRM.STORAGE_DIR, T_PRM.STIM_CSV), ...
    'VariableNamingRule', 'preserve');

window_path = fullfile(T_PRM.OUTPUT_DIR, T_PRM.WINDOW_CSV);
if exist(window_path, 'file')
    windowTbl = readtable(window_path, 'VariableNamingRule', 'preserve');
else
    windowTbl = [];
    fprintf('No window spec CSV (using defaults).\n');
end

if exist(T_PRM.SHANK_MAP_CSV, 'file')
    shankMap = readtable(T_PRM.SHANK_MAP_CSV);
else
    error('ShankMap.csv not found at %s', T_PRM.SHANK_MAP_CSV);
end

i80_path = fullfile(T_PRM.I80_DIR, 'I80_DATA.mat');
if ~exist(i80_path, 'file')
    error('I80_DATA.mat not found at %s', i80_path);
end
load(i80_path, 'I80_DATA');

%% BUILD TIME VECTOR

S = load(fullfile(T_PRM.STORAGE_DIR, T_PRM.RATS{1}, ...
    sprintf('%s_DATA.mat', T_PRM.RATS{1})), 'PRM', 'LF_Shank');
nT_full = size(S.LF_Shank(1).V{1}, 3);
ts_full = -S.PRM.time_extract_preStim  * 1e3;
te_full =  S.PRM.time_extract_postStim * 1e3;
t_full  = linspace(ts_full, te_full, nT_full);
[~, xs] = min(abs(t_full + 50));
[~, xe] = min(abs(t_full - 150));
T_VEC   = t_full(xs:xe);
[~, t_b_start] = min(abs(T_VEC - PARAMS.baseline_window_ms(1)));
[~, t_b_end]   = min(abs(T_VEC - PARAMS.baseline_window_ms(2)));
clear S;

%% RUN DETECTION

SR_DATA = struct();
for iRat = 1:length(T_PRM.RATS)
    rat_id = T_PRM.RATS{iRat};
    age    = T_PRM.AGES{iRat};
    ratFld = sprintf('Rat%s', rat_id);
    fprintf('Rat %s (%s)...\n', rat_id, age);

    fp = fullfile(T_PRM.STORAGE_DIR, rat_id, sprintf('%s_DATA.mat', rat_id));
    if ~exist(fp, 'file'); fprintf('  [SKIP] no _DATA.mat\n'); continue; end

    [L2_shank, L5_shank] = lookup_layer_shanks(rat_id, shankMap);
    if isnan(L2_shank) || isnan(L5_shank)
        fprintf('  [SKIP] L2/L5 not in ShankMap\n'); continue;
    end

    S_rat    = load(fp, 'STIM', 'LF_Shank');
    STIM     = S_rat.STIM;
    LF_Shank = S_rat.LF_Shank;

    SR_DATA.(ratFld).rat_id   = rat_id;
    SR_DATA.(ratFld).age      = age;
    SR_DATA.(ratFld).L2_shank = L2_shank;
    SR_DATA.(ratFld).L5_shank = L5_shank;

    layer_shanks = [L2_shank, L5_shank];
    layer_tags   = {'L2', 'L5'};

    for li = 1:numel(layer_shanks)
        shIdx = layer_shanks(li);
        shFld = sprintf('Shank%d', shIdx);
        SR_DATA.(ratFld).(shFld).layer_designation = layer_tags{li};

        allDepths = LF_Shank(shIdx).SortedDepths;
        IL_chans  = find(allDepths >= PARAMS.PL_IL_boundary_um);
        PL_chans  = find(allDepths <  PARAMS.PL_IL_boundary_um);
        sessions_struct = [];

        for sess = 1:2
            sess_name   = STIM.epochs{sess};
            valid_stims = get_valid_stim_indices(allStimTbl, rat_id, sess_name, ...
                STIM.sorted_current_uA);
            [stim_idx, stim_currents_5, extras_pool] = extract_I80_5stims( ...
                I80_DATA, rat_id, sess_name, valid_stims, STIM.sorted_current_uA);

            sess_entry = struct();
            sess_entry.name            = sess_name;
            sess_entry.stim_idx_sorted = stim_idx(:);
            sess_entry.stim_currents   = stim_currents_5(:);

            Vb_all = LF_Shank(shIdx).V{sess}(:, :, xs:xe);

            [pl_mono_w, pl_poly_w, pl_mono_ok, pl_poly_ok] = ...
                lookup_window_spec_full(windowTbl, rat_id, sess_name, ...
                shIdx, 'PL', PARAMS.PL_mono_window_ms, PARAMS.PL_poly_window_ms);
            sess_entry.PL = run_region_detection_replace(Vb_all, T_VEC, PL_chans, ...
                stim_idx, extras_pool, allDepths, t_b_start, t_b_end, ...
                pl_mono_w, pl_poly_w, pl_mono_ok, pl_poly_ok, PARAMS);

            [il_mono_w, il_poly_w, il_mono_ok, il_poly_ok] = ...
                lookup_window_spec_full(windowTbl, rat_id, sess_name, ...
                shIdx, 'IL', PARAMS.IL_mono_window_ms, PARAMS.IL_poly_window_ms);
            sess_entry.IL = run_region_detection_replace(Vb_all, T_VEC, IL_chans, ...
                stim_idx, extras_pool, allDepths, t_b_start, t_b_end, ...
                il_mono_w, il_poly_w, il_mono_ok, il_poly_ok, PARAMS);

            clear Vb_all;
            if isempty(sessions_struct); sessions_struct = sess_entry;
            else; sessions_struct(end+1, 1) = sess_entry; %#ok<AGROW>
            end
        end
        SR_DATA.(ratFld).(shFld).sessions = sessions_struct;

        %% Fallback for rat 11044 PL on L5
        if strcmp(rat_id, '11044') && shIdx == L5_shank
            fallback_shanks = [3, L2_shank];
            for sess = 1:numel(SR_DATA.(ratFld).(shFld).sessions)
                if SR_DATA.(ratFld).(shFld).sessions(sess).PL.mono.peak_detected
                    continue
                end
                sess_name   = SR_DATA.(ratFld).(shFld).sessions(sess).name;
                valid_stims = get_valid_stim_indices(allStimTbl, rat_id, sess_name, ...
                    STIM.sorted_current_uA);
                [stim_idx, ~, extras_pool] = extract_I80_5stims(I80_DATA, rat_id, sess_name, ...
                    valid_stims, STIM.sorted_current_uA);
                for fb = 1:numel(fallback_shanks)
                    fb_sh = fallback_shanks(fb);
                    if fb_sh == shIdx; continue; end
                    fb_depths = LF_Shank(fb_sh).SortedDepths;
                    fb_PL     = find(fb_depths < PARAMS.PL_IL_boundary_um);
                    Vb_fb     = LF_Shank(fb_sh).V{sess}(:, :, xs:xe);
                    [pl_mono_w, pl_poly_w, pl_mono_ok, pl_poly_ok] = ...
                        lookup_window_spec_full(windowTbl, rat_id, sess_name, ...
                        fb_sh, 'PL', PARAMS.PL_mono_window_ms, PARAMS.PL_poly_window_ms);
                    fb_PL_det = run_region_detection_replace(Vb_fb, T_VEC, fb_PL, ...
                        stim_idx, extras_pool, fb_depths, t_b_start, t_b_end, ...
                        pl_mono_w, pl_poly_w, pl_mono_ok, pl_poly_ok, PARAMS);
                    clear Vb_fb;
                    if fb_PL_det.mono.peak_detected
                        SR_DATA.(ratFld).(shFld).sessions(sess).PL = fb_PL_det;
                        break
                    end
                end
            end
        end
    end
    clear S_rat LF_Shank STIM;
end

save(fullfile(T_PRM.OUTPUT_DIR, 'SR_DATA_peaks.mat'), ...
    'SR_DATA', 'T_VEC', 'PARAMS', 'T_PRM', '-v7.3');
fprintf('\nSaved: SR_DATA_peaks.mat\n');

i80_table_path = fullfile(T_PRM.I80_DIR, 'I80_TABLE_IL_Mono_CLEAN.csv');
export_I80_table_clean(I80_DATA, i80_table_path);

%% VERIFICATION PLOTS

plot_slope_verification(SR_DATA, T_VEC, PARAMS, T_PRM.FIG_DIR);
fprintf('\nMANUALLY CHECK SLOPE DETECTION:\n  %s\n', T_PRM.FIG_DIR);

%% AVERAGE-TRACE FIGURES

plot_avg_traces_by_age(SR_DATA, T_VEC, T_PRM, PARAMS, ...
    'PRM_amplitude_var', 'lowest');

%% PEAK WINDOW TIMES CSV

peak_window_csv = fullfile(T_PRM.OUTPUT_DIR, 'peak_window_times_post_i80.csv');
build_peak_window_times_csv(SR_DATA, T_VEC, T_PRM, peak_window_csv);

%% CHANNEL INCLUSION VERIFICATION

chan_depth_dir = fullfile(T_PRM.OUTPUT_DIR, 'ChannelDepth_Violins');
plot_channel_depth_violins(SR_DATA, chan_depth_dir, 'PRM_threshold', '50pct');
plot_channel_depth_violins(SR_DATA, chan_depth_dir, 'PRM_threshold', '25pct');

%% PEAK SUMMARY TABLES

peak_csv = fullfile(T_PRM.OUTPUT_DIR, 'PeakSummary_Mono_35ms.csv');
build_peak_summary_table(SR_DATA, ...
    'PRM_windows',  {'mono'}, ...
    'PRM_regions',  {'PL','IL'}, ...
    'PRM_save_csv', peak_csv);

build_peak_summary_table(SR_DATA, ...
    'PRM_windows',  {'mono','poly'}, ...
    'PRM_save_csv', fullfile(T_PRM.OUTPUT_DIR, 'PeakSummary_All_35ms.csv'));

fprintf('Peak summary CSV: %s\n', peak_csv);
