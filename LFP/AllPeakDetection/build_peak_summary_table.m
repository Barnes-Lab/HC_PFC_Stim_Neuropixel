function T = build_peak_summary_table(SR_DATA, varargin)

% Per-rat peak summary table Rows with no detection ie PL monoskipped


%% PARAMETERS TO CHANGE
PRM_windows  = {'mono'};
PRM_regions  = {'PL', 'IL'};
PRM_save_csv = '';
Extract_varargin;

layers = {'L2', 'L5'};
rats   = fieldnames(SR_DATA);
rats   = rats(startsWith(rats, 'Rat'));

n_max  = numel(rats) * numel(layers) * 2 * numel(PRM_regions) * numel(PRM_windows);
rows   = cell(n_max, 1);
n_used = 0;

for iR = 1:numel(rats)
    rd = SR_DATA.(rats{iR});
    if ~isfield(rd, 'L2_shank'); continue; end
    rat_str = string(rd.rat_id);
    age_str = string(rd.age);
    layer_shanks = [rd.L2_shank, rd.L5_shank];

    for li = 1:numel(layers)
        shFld = sprintf('Shank%d', layer_shanks(li));
        if ~isfield(rd, shFld); continue; end
        sess_struct = rd.(shFld).sessions;

        for sess = 1:numel(sess_struct)
            sess_name = string(sess_struct(sess).name);
            for ir = 1:numel(PRM_regions)
                reg = PRM_regions{ir};
                if ~isfield(sess_struct(sess), reg); continue; end
                rdata = sess_struct(sess).(reg);
                for iw = 1:numel(PRM_windows)
                    wk = PRM_windows{iw};
                    if ~isfield(rdata, wk); continue; end
                    pk = rdata.(wk);
                    if ~isstruct(pk) || ~isfield(pk, 'peak_detected'); continue; end
                    if ~pk.peak_detected; continue; end

                    n_used = n_used + 1;
                    rows{n_used} = make_row(rat_str, age_str, sess_name, ...
                        string(layers{li}), string(reg), string(wk), pk);
                end
            end
        end
    end
end

if n_used == 0
    T = table(); 
    return;
end

T = vertcat(rows{1:n_used});

if ~isempty(PRM_save_csv)
    [out_dir, ~, ~] = fileparts(PRM_save_csv);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir'); mkdir(out_dir); end
    writetable(T, PRM_save_csv);
end

end

% Easier to have in line
function R = make_row(rat, age, sess, layer, reg, win, pk)

amp_anchor   = mean(pk.best_channel_amplitude, 'omitnan');
lat_anchor   = mean(pk.best_channel_latency,   'omitnan');
amp_50       = mean(pk.mean_amplitude_50pct,   'omitnan');
lat_50       = mean(pk.mean_latency_50pct,     'omitnan');
amp_25       = mean(pk.mean_amplitude_25pct,   'omitnan');
lat_25       = mean(pk.mean_latency_25pct,     'omitnan');
n_responding = sum(~isnan(pk.best_channel_amplitude));

R = table( rat, age, sess, layer, reg, win, ...
    pk.best_channel, pk.best_channel_depth, ...
    pk.template_amplitude, pk.template_latency, ...
    amp_anchor, lat_anchor, ...
    amp_50, lat_50, ...
    amp_25, lat_25, ...
    pk.n_region_channels_50pct, pk.n_region_channels_25pct, ...
    n_responding, ...
    'VariableNames', {'Rat','Age','Session','Layer','Region','Window', ...
                      'Best_Channel','Best_Channel_Depth_um', ...
                      'Template_Amp_mV','Template_Lat_ms', ...
                      'Anchor_Amp_mean_mV','Anchor_Lat_mean_ms', ...
                      'Mean_Amp_50pct_mV','Mean_Lat_50pct_ms', ...
                      'Mean_Amp_25pct_mV','Mean_Lat_25pct_ms', ...
                      'N_chans_50pct_template','N_chans_25pct_template', ...
                      'N_stims_responding'});

end
