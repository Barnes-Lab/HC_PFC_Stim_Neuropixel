function [data, stim_report] = extract_negative_data( ...
    SR_DATA, T_PRM, amplitude_var, layer_filter)

%
% Extract per-(rat, layer, session, region, window) negative-peak
% amplitude and latency from SR_DATA, 
%  The 5 analysis stims per recording are
% pre-selected upstream by the I80 detection stage; this function
% just averages whatever per-stim arrays detect_neg_peak_5stim stored.
%
% SS 2026

%%  INPUT PARAMS

variant = amplitude_var;
[amp_field, lat_field] = variant_fields(variant);

fprintf('=== Extracting negative-peak data ===\n');
fprintf('Variant      : %s  (amp=%s, lat=%s)\n', variant, amp_field, lat_field);
fprintf('Layer filter : %s\n', layer_filter);

%%  LAYOUT
if strcmp(layer_filter, 'Combined')
    layers_master = {'L2', 'L5'};
else
    layers_master = {layer_filter};
end
n_layers_per_rat = numel(layers_master);
n_entries        = numel(T_PRM.RATS) * n_layers_per_rat;

data           = struct();
data.amplitude = nan(n_entries, 2, 2, 2);
data.latency   = nan(n_entries, 2, 2, 2);
data.rat_ids   = cell(n_entries, 1);
data.ages      = cell(n_entries, 1);
data.layers    = cell(n_entries, 1);
data.young_idx = [];
data.old_idx   = [];

report_entries = {};
entry_idx = 0;

%%  LOOP rats x layers x sessions x regions x windows
for iRat = 1:numel(T_PRM.RATS)
    rat      = T_PRM.RATS{iRat};
    age      = T_PRM.AGES{iRat};
    ratField = sprintf('Rat%s', rat);

    if ~isfield(SR_DATA, ratField)
        fprintf('  [WARN] %s not found in SR_DATA\n', rat);
        continue;
    end

    L2_shank = SR_DATA.(ratField).L2_shank;
    L5_shank = SR_DATA.(ratField).L5_shank;

    for iLay = 1:n_layers_per_rat
        layer     = layers_master{iLay};
        entry_idx = entry_idx + 1;

        data.rat_ids{entry_idx} = rat;
        data.ages{entry_idx}    = age;
        data.layers{entry_idx}  = layer;
        if strcmp(age, 'Y')
            data.young_idx(end+1) = entry_idx; %#ok<AGROW>
        else
            data.old_idx(end+1)   = entry_idx; %#ok<AGROW>
        end

        if strcmp(layer, 'L2'); shank_num = L2_shank;
        else;                   shank_num = L5_shank;
        end
        shLabel = sprintf('Shank%d', shank_num);
        if ~isfield(SR_DATA.(ratField), shLabel)
            fprintf('  [WARN] %s %s not found\n', rat, shLabel);
            continue;
        end

        sessions_struct = SR_DATA.(ratField).(shLabel).sessions;
        n_sess = min(2, numel(sessions_struct));

        for sess = 1:n_sess
            session_data = sessions_struct(sess);
            if contains(session_data.name, 'iHC', 'IgnoreCase', true)
                stim_site = 'iHC';
            else
                stim_site = 'vHC';
            end

            for r = 1:2
                if r == 1; reg_name = 'PL'; else; reg_name = 'IL'; end
                if ~isfield(session_data, reg_name); continue; end

                report_entries(end+1, :) = {rat, age, shank_num, layer, ...
                    session_data.name, stim_site, reg_name, variant}; %#ok<AGROW>

                for w = 1:2
                    if w == 1; win_name = 'mono'; else; win_name = 'poly'; end
                    if ~isfield(session_data.(reg_name), win_name); continue; end
                    win_data = session_data.(reg_name).(win_name);

                    [amp, lat] = extract_mean_for_variant( ...
                        win_data, amp_field, lat_field);
                    data.amplitude(entry_idx, w, r, sess) = amp;
                    data.latency  (entry_idx, w, r, sess) = lat;
                end
            end
        end
    end
end

%% TRIM 
data.amplitude = data.amplitude(1:entry_idx, :, :, :);
data.latency   = data.latency  (1:entry_idx, :, :, :);
data.rat_ids   = data.rat_ids(1:entry_idx);
data.ages      = data.ages   (1:entry_idx);
data.layers    = data.layers (1:entry_idx);
data.young_idx = data.young_idx(data.young_idx <= entry_idx);
data.old_idx   = data.old_idx (data.old_idx    <= entry_idx);

%% PL NO-RESPONSE set to 0 mV amplitude 
PL_amp = data.amplitude(:, :, 1, :);
PL_amp(isnan(PL_amp)) = 0;
data.amplitude(:, :, 1, :) = PL_amp;

%%  REPORT
stim_report = cell2table(report_entries, 'VariableNames', ...
    {'Rat', 'Age', 'Shank', 'Layer', 'Session', 'Stim_Site', ...
     'Region', 'Variant'});

fprintf('Extracted %d entries (Young=%d, Old=%d) using variant ''%s''\n', ...
    entry_idx, numel(data.young_idx), numel(data.old_idx), variant);

end




