function create_depth_psth_from_comprehensive(all_SP_filtered, ALL_STIM_TIMES, sess_types, ageCols, figDir)
% Depth-resolved PSTH per session
% type plus a combined-sessions panel.
%  For each session in
% [sess_types, 'Combined'], averages per-rat depth PSTHs across each age
+

%
% SS 2026 adapted from steinmetz

depthBinSize = 30;
timeBinSize  = 5;
winResp      = [-50 150];
bslWin       = [-400 0];
PL_depth     = 3600;

largest_bin      = 5125;
smallest_bin     = largest_bin - (96 * depthBinSize);
global_depthBins = smallest_bin:depthBinSize:largest_bin;

sess_types_all = [sess_types, {'Combined'}];

for s_idx = 1:length(sess_types_all)
    s_type = sess_types_all{s_idx};
    fprintf('  Processing %s...\n', s_type);

    Y_sum = []; O_sum = []; Y_count = 0; O_count = 0;
    unique_rats = unique({all_SP_filtered.rat});

    for r_idx = 1:length(unique_rats)
        rat_id  = unique_rats{r_idx};
        rat_SP  = all_SP_filtered(strcmp({all_SP_filtered.rat}, rat_id));
        rat_age = rat_SP(1).age;

        event_times = collect_rat_events(ALL_STIM_TIMES, rat_id, s_type, sess_types);
        if isempty(event_times); continue; end

        [spikeTimes, spikeDepths] = pool_rat_spikes(rat_SP);
        [timeBins, ~, allP_norm]  = psthByDepthNormalized(spikeTimes, spikeDepths, ...
            global_depthBins, timeBinSize/1000, event_times, winResp/1000, bslWin/1000);
        allP_norm = flipud(allP_norm);

        if strcmp(rat_age, 'Y')
            if isempty(Y_sum); Y_sum = allP_norm; else; Y_sum = Y_sum + allP_norm; end
            Y_count = Y_count + 1;
        else
            if isempty(O_sum); O_sum = allP_norm; else; O_sum = O_sum + allP_norm; end
            O_count = O_count + 1;
        end
    end
    if Y_count == 0 || O_count == 0; continue; end

    allP_Y      = Y_sum / Y_count;
    allP_O      = O_sum / O_count;
    timeBins_ms = timeBins * 1000;

    % Heatmap: Y vs O with shared color scale
    figure('Position', [100, 100, 1600, 600]);
    subplot(1, 2, 1);
    plot_psth_depth_heatmap(timeBins_ms, global_depthBins, allP_Y, ...
        sprintf('Young (n=%d)', Y_count), PL_depth);
    subplot(1, 2, 2);
    plot_psth_depth_heatmap(timeBins_ms, global_depthBins, allP_O, ...
        sprintf('Old (n=%d)', O_count), PL_depth);

    max_val = max([allP_Y(:); allP_O(:)]);
    if     max_val < 10; clim_max = ceil(max_val);
    elseif max_val < 50; clim_max = ceil(max_val/5) * 5;
    else;                clim_max = ceil(max_val/10) * 10;
    end
    subplot(1, 2, 1); clim([0, clim_max]);
    subplot(1, 2, 2); clim([0, clim_max]);

    sgtitle(sprintf('Depth-Resolved PSTH: %s', s_type), 'FontSize', 16);
    saveas(gcf, fullfile(figDir, sprintf('depth_psth_%s.png', s_type)));
    saveas(gcf, fullfile(figDir, sprintf('depth_psth_%s.fig', s_type)));
    close(gcf);

    % PL vs IL traces per age, shared y-axis
    depth_centers = global_depthBins(1:end-1);
    PL_idx = depth_centers <  PL_depth;
    IL_idx = depth_centers >= PL_depth;

    figure('Position', [100, 100, 1400, 600]);
    subplot(2, 2, 1); plot_region_psth_trace(timeBins_ms, allP_Y(PL_idx, :), 'Young PL', ageCols.Y);
    subplot(2, 2, 3); plot_region_psth_trace(timeBins_ms, allP_Y(IL_idx, :), 'Young IL', ageCols.Y);
    subplot(2, 2, 2); plot_region_psth_trace(timeBins_ms, allP_O(PL_idx, :), 'Old PL',   ageCols.O);
    subplot(2, 2, 4); plot_region_psth_trace(timeBins_ms, allP_O(IL_idx, :), 'Old IL',   ageCols.O);

    y_max = 0;
    for i = 1:4; subplot(2, 2, i); y_max = max(y_max, max(ylim)); end
    for i = 1:4; subplot(2, 2, i); ylim([0, y_max]); end

    sgtitle(sprintf('PL vs IL Response: %s', s_type), 'FontSize', 16);
    saveas(gcf, fullfile(figDir, sprintf('depth_psth_PL_IL_%s.png', s_type)));
    saveas(gcf, fullfile(figDir, sprintf('depth_psth_PL_IL_%s.fig', s_type)));
    close(gcf);
end
fprintf('  Depth PSTH plots created!\n');
end


function event_times = collect_rat_events(ALL_STIM_TIMES, rat_id, s_type, sess_types)
% Pull event times for one rat. s_type='Combined' pools across sess_types.
if strcmp(s_type, 'Combined'); sess_list = sess_types;
else;                          sess_list = {s_type};
end
event_times = [];
for s = 1:length(sess_list)
    T = ALL_STIM_TIMES(strcmp(string(ALL_STIM_TIMES.Rat), rat_id) & ...
                       strcmp(ALL_STIM_TIMES.Session, sess_list{s}), :);
    if ismember('keep', T.Properties.VariableNames); T = T(~strcmp(T.keep, 'N'), :); end
    if ~isempty(T); event_times = [event_times; T.Time_s]; end
end
end


function [spikeTimes, spikeDepths] = pool_rat_spikes(rat_SP)
% Pool good spikes and per-spike depths across one rat's units.
spikeTimes  = [];
spikeDepths = [];
for u = 1:numel(rat_SP)
    if isfield(rat_SP(u), 'IX_of_good_spikes')
        good_idx = logical(rat_SP(u).IX_of_good_spikes);
        if length(good_idx) > length(rat_SP(u).t_uS)
            good_idx = good_idx(1:length(rat_SP(u).t_uS));
        elseif length(good_idx) < length(rat_SP(u).t_uS)
            tmp = false(size(rat_SP(u).t_uS));
            tmp(1:length(good_idx)) = good_idx;
            good_idx = tmp;
        end
    else
        good_idx = true(size(rat_SP(u).t_uS));
    end
    unit_spike_times = rat_SP(u).t_uS(good_idx) / 1e6;
    unit_depths      = repmat(rat_SP(u).neuropixels_depth_uM, sum(good_idx), 1);
    spikeTimes       = [spikeTimes; unit_spike_times(:)];
    spikeDepths      = [spikeDepths; unit_depths(:)];
end
end
