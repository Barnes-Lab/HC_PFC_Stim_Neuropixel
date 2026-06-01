function plot_per_rat_pyr_depth_psth(all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
    I80_DATA, sess_types, PL_depth, rats, figDir, varargin)
% PLOT_PER_RAT_PYR_DEPTH_PSTH  Per-rat sanity check. For each rat: 1 PNG with
% 4 panels (iHC, vHC) x (PL, IL). Each panel = depth-sorted heatmap of
% baseline-subtracted PSTHs across PYR neurons, time-locked to the 5
% I80-selected stims. Rows = PYR neurons sorted by depth. Color = dFR (Hz).
% Color limits shared across the 4 panels of each rat.
% SS 2026

PRM_TIME_WIN_S = [-0.1, 0.15];
PRM_BSL_WIN_S  = [-0.1, -0.01];
PRM_BIN_SIZE_S = 0.005;
Extract_varargin;

regions     = {'PL', 'IL'};
time_bins   = PRM_TIME_WIN_S(1) : PRM_BIN_SIZE_S : PRM_TIME_WIN_S(2);
bin_centers = time_bins(1:end-1) + PRM_BIN_SIZE_S/2;
bsl_idx     = (bin_centers >= PRM_BSL_WIN_S(1)) & (bin_centers < PRM_BSL_WIN_S(2));

if ~exist(figDir, 'dir'); mkdir(figDir); end
cm = build_diverging_cmap();

for r = 1:numel(rats)
    rat = rats{r};
    psth_cell = collect_rat_psth(all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
        I80_DATA, rat, sess_types, regions, PL_depth, time_bins, bsl_idx);
    if isempty(fieldnames(psth_cell)); continue; end

    all_v = [];
    for s = 1:numel(sess_types)
        for rg = 1:numel(regions)
            M = psth_cell.(sess_types{s}).(regions{rg}).M;
            if ~isempty(M); all_v = [all_v; M(:)]; end
        end
    end
    if isempty(all_v); continue; end
    cl_max = max(abs(all_v), [], 'omitnan');
    if isempty(cl_max) || cl_max == 0 || isnan(cl_max); cl_max = 1; end
    cl = [-cl_max, cl_max];

    figure('Position', [100, 100, 1100, 900]);
    tiledlayout(numel(sess_types), numel(regions), 'TileSpacing', 'compact', 'Padding', 'compact');
    for s = 1:numel(sess_types)
        for rg = 1:numel(regions)
            ax  = nexttile;
            p   = psth_cell.(sess_types{s}).(regions{rg});
            ttl = sprintf('%s - %s (n=%d)', sess_types{s}, regions{rg}, size(p.M, 1));
            if isempty(p.M)
                text(ax, 0.5, 0.5, 'no data', 'HorizontalAlignment', 'center', 'Units', 'normalized');
                axis(ax, 'off'); title(ax, ttl); continue;
            end
            imagesc(ax, bin_centers * 1000, 1:size(p.M, 1), p.M);
            set(ax, 'YDir', 'reverse');
            clim(ax, cl);
            colormap(ax, cm);
            hold(ax, 'on');
            xline(ax, 0, 'k-', 'LineWidth', 1.5);
            xlabel(ax, 'Time from stim (ms)');
            ylabel(ax, 'PYR neurons (depth-sorted)');
            title(ax, ttl);
            cb = colorbar(ax);
            cb.Label.String = '\Delta FR (Hz)';
        end
    end
    sgtitle(sprintf('Rat %s - PYR \\Delta FR depth-PSTH (5 I80 stims)', rat));
    saveas(gcf, fullfile(figDir, sprintf('depth_psth_PYR_%s.png', rat)));
    close(gcf);
end

fprintf('  Per-rat PYR depth-PSTH: %s\n', figDir);
end


function out = collect_rat_psth(all_SP, T_stims, classif, I80_DATA, ...
    rat, sess_types, regions, PL_depth, time_bins, bsl_idx)
out = struct();
rats_col = {all_SP.rat}';
n_total  = numel(all_SP);
is_keep  = false(n_total, 1);
for u = 1:n_total
    is_keep(u) = strcmp(classif.type{u}, 'PYR') && strcmp(rats_col{u}, rat);
end
ix = find(is_keep);
if isempty(ix); return; end

n = numel(ix);
neur_depth  = nan(n, 1);
neur_region = cell(n, 1);
for k = 1:n
    SP = all_SP(ix(k));
    if isfield(SP, 'neuropixels_depth_uM');     neur_depth(k) = SP.neuropixels_depth_uM;
    elseif isfield(SP, 'depth_uM');             neur_depth(k) = SP.depth_uM;
    end
    if isnan(neur_depth(k)); continue; end
    if neur_depth(k) < PL_depth; neur_region{k} = 'PL';
    else;                        neur_region{k} = 'IL';
    end
end

for s = 1:numel(sess_types)
    sess = sess_types{s};
    ev   = get_5_I80_events(T_stims, I80_DATA, rat, sess);
    for rg = 1:numel(regions)
        region = regions{rg};
        sel = strcmp(neur_region, region);
        if isempty(ev) || ~any(sel)
            out.(sess).(region).M      = [];
            out.(sess).(region).depths = [];
            continue;
        end
        ix_sel = ix(sel);
        D      = neur_depth(sel);
        nB     = numel(time_bins) - 1;
        M      = nan(numel(ix_sel), nB);
        for k = 1:numel(ix_sel)
            SP   = all_SP(ix_sel(k));
            spk  = double(SP.t_uS) / 1e6;
            psth = calculate_psth(spk, ev, time_bins);
            bsl  = mean(psth(bsl_idx), 'omitnan');
            M(k, :) = psth - bsl;
        end
        [D_sorted, sord] = sort(D);
        out.(sess).(region).M      = M(sord, :);
        out.(sess).(region).depths = D_sorted;
    end
end
end


function ev = get_5_I80_events(T_stims, I80_DATA, rat, sess)
T = T_stims(strcmp(string(T_stims.Rat), rat) & strcmp(T_stims.Session, sess), :);
if isempty(T); ev = []; return; end
cur_col = '';
for cc = {'Current_uA', 'Intensity_uA', 'Stim_uA'}
    if ismember(cc{1}, T.Properties.VariableNames); cur_col = cc{1}; break; end
end
if isempty(cur_col); ev = []; return; end
[ev, ~, ~] = select_I80_event_times_neuron(I80_DATA, rat, sess, T, cur_col);
end


function cm = build_diverging_cmap()
n   = 128;
neg = [linspace(0.13, 1, n)', linspace(0.13, 1, n)', linspace(0.55, 1, n)'];
pos = [linspace(1, 0.70, n)', linspace(1, 0,    n)', linspace(1, 0,    n)'];
cm  = [neg; pos];
end
