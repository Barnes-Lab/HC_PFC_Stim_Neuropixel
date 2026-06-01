%% ========================================================================
% RUN_LF_MONO_DETECTION
% ========================================================================
% Mono LFP peak detection in IL. D1 (lowest template min) then
% V2 (per-channel findpeaks) 

clear; close all; clc;

%% INPUT PARAMETERS


T_PRM.RATS={'10940','10947','10986','10987','10994',...
    '11039','11042','11043','11040','11044'};
T_PRM.AGES={'O','O','Y','Y','O','O','Y','Y','O','Y'};

T_PRM.STORAGE_DIR   = 'C:\DATA\NP2\Processed';
T_PRM.OUTPUT_DIR    = fullfile(T_PRM.STORAGE_DIR, 'LF_Stimulus_Response_Analysis');
T_PRM.FIG_DIR       = fullfile(T_PRM.OUTPUT_DIR, 'MonoDetection_Verify');
T_PRM.STIM_CSV      = 'ALL_STIM_TIMES_KEEP_CHECK.csv';
T_PRM.WINDOW_CSV    = 'manual_window_specifications_mono.csv';
T_PRM.SHANK_MAP_CSV = fullfile(T_PRM.STORAGE_DIR, 'ShankInfo', 'ShankMap.csv');

PARAMS = struct();
PARAMS.IL_mono_window_ms     = [10, 40];
PARAMS.PL_IL_boundary_um     = 3600;
PARAMS.baseline_window_ms    = [-410, -10];
PARAMS.prepeak_baseline_ms   = [5, 9];
PARAMS.smooth_window_ms      = 3;
PARAMS.drift_thresh_mV       = 0.15;
PARAMS.n_strongest           = 3;
PARAMS.wall_tol_ms           = 1.0;
PARAMS.wall_hit_frac         = 0.5;
PARAMS.prominence_sigma      = 3;
PARAMS.min_peak_distance_ms  = 5;
PARAMS.min_peak_width_ms     = 0.5;
PARAMS.per_trial_half_win_ms = 10;
PARAMS.colors.IL             = [0.88, 0.35, 0.11];

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
    fprintf('No window spec CSV (using default %g-%g ms).\n', PARAMS.IL_mono_window_ms);
end

if exist(T_PRM.SHANK_MAP_CSV, 'file')
    shankMap = readtable(T_PRM.SHANK_MAP_CSV);
else
    shankMap = [];
    fprintf('WARNING ShankMap.csv not found no layer specifications.\n');
end

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

MONO_DATA = struct();
n_d1 = 0; n_v2 = 0; n_fail_v2 = 0;

for iRat = 1:length(T_PRM.RATS)
    rat_id = T_PRM.RATS{iRat};
    age    = T_PRM.AGES{iRat};
    ratFld = sprintf('Rat%s', rat_id);
    fprintf('Rat %s (%s)...\n', rat_id, age);

    fp = fullfile(T_PRM.STORAGE_DIR, rat_id, sprintf('%s_DATA.mat', rat_id));
    if ~exist(fp, 'file'); fprintf('  [SKIP] no _DATA.mat\n'); continue; end

    S_rat    = load(fp, 'STIM', 'LF_Shank', 'PRM');
    STIM     = S_rat.STIM;
    LF_Shank = S_rat.LF_Shank;
    PRM_rat  = S_rat.PRM;

    MONO_DATA.(ratFld).rat_id = rat_id;
    MONO_DATA.(ratFld).age    = age;
    [L2_shank, L5_shank] = lookup_layer_shanks(rat_id, shankMap);
    MONO_DATA.(ratFld).L2_shank = L2_shank;
    MONO_DATA.(ratFld).L5_shank = L5_shank;

    for shIdx = 1:PRM_rat.shank_no
        shFld = sprintf('Shank%d', shIdx);
        if     shIdx == L2_shank; layer_tag = 'L2';
        elseif shIdx == L5_shank; layer_tag = 'L5';
        else;                     layer_tag = '';
        end
        MONO_DATA.(ratFld).(shFld).layer_designation = layer_tag;

        allDepths = LF_Shank(shIdx).SortedDepths;
        IL_chans  = find(allDepths >= PARAMS.PL_IL_boundary_um);
        sessions_struct = [];

        for sess = 1:2
            sess_name   = STIM.epochs{sess};
            valid_stims = get_valid_stim_indices(allStimTbl, rat_id, sess_name, ...
                STIM.sorted_current_uA);
            if isempty(valid_stims); continue; end

            [mono_win, mono_valid, has_ov] = lookup_window_spec( ...
                windowTbl, rat_id, sess_name, shIdx, 'IL', ...
                PARAMS.IL_mono_window_ms);

            Vb_all = LF_Shank(shIdx).V{sess}(:, :, xs:xe);

            sess_entry = struct();
            sess_entry.name            = sess_name;
            sess_entry.stim_currents   = STIM.sorted_current_uA(:);
            sess_entry.valid_stims     = valid_stims(:);
            sess_entry.mono_window     = mono_win;
            sess_entry.window_override = has_ov;
            sess_entry.IL = init_il_entry(numel(valid_stims), numel(IL_chans), mono_win);

            if isempty(IL_chans) || ~mono_valid
                sess_entry.IL.channel_mode = 'invalid';
            else
                pre = preprocess_lfp_traces(Vb_all, T_VEC, IL_chans, valid_stims, ...
                    t_b_start, t_b_end, ...
                    'PRM_prepeak_baseline_ms', PARAMS.prepeak_baseline_ms, ...
                    'PRM_smooth_window_ms',    PARAMS.smooth_window_ms, ...
                    'PRM_drift_thresh_mV',     PARAMS.drift_thresh_mV, ...
                    'PRM_strong_window_ms',    mono_win, ...
                    'PRM_n_strongest',         PARAMS.n_strongest);

                pk_d1 = detect_mono_d1(pre.V_proc, T_VEC, IL_chans, valid_stims, ...
                    pre.candidate_local, allDepths, mono_win, ...
                    'PRM_n_strongest', PARAMS.n_strongest);

                [is_walled, frac_walled] = check_window_wall_hit( ...
                    pk_d1, mono_win, PARAMS.wall_tol_ms, PARAMS.wall_hit_frac);

                final_pk = pk_d1;
                mode_str = 'D1';

                if is_walled
                    pk_v2 = detect_mono_v2_findpeaks(pre.V_proc, T_VEC, IL_chans, ...
                        valid_stims, pre.candidate_local, allDepths, mono_win, ...
                        t_b_start, t_b_end, ...
                        'PRM_n_strongest',           PARAMS.n_strongest, ...
                        'PRM_prominence_sigma',      PARAMS.prominence_sigma, ...
                        'PRM_min_peak_distance_ms',  PARAMS.min_peak_distance_ms, ...
                        'PRM_min_peak_width_ms',     PARAMS.min_peak_width_ms, ...
                        'PRM_per_trial_half_win_ms', PARAMS.per_trial_half_win_ms);
                    if pk_v2.peak_detected
                        final_pk = pk_v2;
                        mode_str = 'V2';
                        n_v2 = n_v2 + 1;
                    else
                        mode_str  = 'V2_failed_keep_D1';
                        n_fail_v2 = n_fail_v2 + 1;
                    end
                else
                    n_d1 = n_d1 + 1;
                end

                sess_entry.IL.mono              = final_pk;
                sess_entry.IL.strong_stims      = pre.strong_stims;
                sess_entry.IL.excluded_channels = pre.excluded_global;
                sess_entry.IL.channel_mode      = mode_str;
                sess_entry.IL.frac_walled       = frac_walled;
            end

            if isempty(sessions_struct); sessions_struct = sess_entry;
            else; sessions_struct(end+1, 1) = sess_entry; %#ok<AGROW>
            end
            clear Vb_all sess_entry pre;
        end
        MONO_DATA.(ratFld).(shFld).sessions = sessions_struct;
    end
    clear S_rat LF_Shank STIM;
end

fprintf('\nDetection summary: D1=%d  V2_fallback=%d  V2_failed=%d\n', ...
    n_d1, n_v2, n_fail_v2);

save(fullfile(T_PRM.OUTPUT_DIR, 'MONO_DATA.mat'), ...
    'MONO_DATA', 'T_VEC', 'PARAMS', '-v7.3');
fprintf('Saved: MONO_DATA.mat\n');

%% VERIFICATION PLOTS

plot_mono_verification(MONO_DATA, T_VEC, PARAMS, T_PRM.FIG_DIR);
fprintf('\nMANUALLY CHECK MONO DETECTION:\n  %s\n', T_PRM.FIG_DIR);
