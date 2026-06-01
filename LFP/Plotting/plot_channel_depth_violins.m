function plot_channel_depth_violins(SR_DATA, FIG_DIR, varargin)
% Channels with a response shown 
%
%  anatomical region range bounds using ksdensity

%
% SS 2026

PRM_threshold = '50pct';
PRM_color_Y   = [0.20 0.55 0.20];
PRM_color_O   = [0.50 0.20 0.55];
Extract_varargin;

if ~ismember(PRM_threshold, {'50pct', '25pct'})
    error('PRM_threshold must be ''50pct'' or ''25pct''.');
end
if ~exist(FIG_DIR, 'dir'); mkdir(FIG_DIR); end

T = build_channel_inclusion_long(SR_DATA, PRM_threshold);
if isempty(T); fprintf('No channel inclusion data to plot.\n'); return; end

depth_bounds = compute_region_depth_bounds(SR_DATA);

sessions = {'iHC', 'vHC'};
windows  = {'mono', 'poly'};
win_lbls = {'Mono', 'Poly'};
regions  = {'PL', 'IL'};
layers   = {'L2', 'L5'};

for is = 1:2
    for iw = 1:2
        sn = sessions{is};
        wk = windows{iw};
        fig = figure('Position', [50 50 1100 800], 'Color', 'w', 'Visible', 'off');
        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        for ir = 1:2
            for il = 1:2
                nexttile;
                plot_one_violin_panel(T, layers{il}, regions{ir}, wk, sn, ...
                    PRM_color_Y, PRM_color_O, depth_bounds.(regions{ir}));
                title(sprintf('%s - %s', regions{ir}, layers{il}), 'FontSize', 10);
            end
        end
        sgtitle(sprintf('Channel depths (%s) - %s %s', PRM_threshold, sn, win_lbls{iw}), ...
            'FontSize', 12);
        saveas(fig, fullfile(FIG_DIR, sprintf('ChanDepth_%s_%s_%s.png', sn, win_lbls{iw}, PRM_threshold)));
        close(fig);
    end
end
fprintf('  Channel-depth violins (%s) saved to: %s\n', PRM_threshold, FIG_DIR);
end


function plot_one_violin_panel(T, layer, reg, win_kind, sess, color_Y, color_O, region_bounds)
hold on;
sub = T(T.Layer == layer & T.Region == reg & T.Window == win_kind & T.Session == sess, :);
draw_violin(sub.Depth_um(sub.Age == "Y"), 1, color_Y, region_bounds);
draw_violin(sub.Depth_um(sub.Age == "O"), 2, color_O, region_bounds);
set(gca, 'YDir', 'reverse', 'XTick', [1 2], 'XTickLabel', {'Y', 'O'}, 'FontSize', 9);
xlim([0.4 2.6]); ylim(region_bounds);
ylabel('Depth (\mum)', 'FontSize', 9); box off;
end


function draw_violin(depths, x_center, col, bounds)
% Kernel-density violin centred at x_center, bounded to anatomical range
% via ksdensity reflection boundary correction. Mean +/- SEM marked black.
if isempty(depths) || all(isnan(depths)); return; end
d = depths(~isnan(depths));
if numel(d) < 2
    plot(x_center, d, 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 6); return;
end

eps_in  = 1e-6 * (bounds(2) - bounds(1));
support = [bounds(1) - eps_in, bounds(2) + eps_in];
[f, xi] = ksdensity(d, 'Support', support, 'BoundaryCorrection', 'reflection');
f = f / max(f) * 0.35;

fill([x_center - f, fliplr(x_center + f)], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.45, 'EdgeColor', col, 'LineWidth', 1.2);

m  = mean(d);
se = std(d) / sqrt(numel(d));
plot([x_center - 0.18, x_center + 0.18], [m, m], 'k-', 'LineWidth', 2);
plot([x_center, x_center], [m - se, m + se], 'k-', 'LineWidth', 1.5);
end


function bounds = compute_region_depth_bounds(SR_DATA)
% Cohort-wide anatomical range per region: min/max of region_chan_depths_um.
regions = {'PL', 'IL'};
bounds  = struct('PL', [Inf, -Inf], 'IL', [Inf, -Inf]);

rats = fieldnames(SR_DATA);
rats = rats(startsWith(rats, 'Rat'));
for iR = 1:numel(rats)
    rd = SR_DATA.(rats{iR});
    shanks = fieldnames(rd);
    shanks = shanks(startsWith(shanks, 'Shank'));
    for is = 1:numel(shanks)
        if ~isfield(rd.(shanks{is}), 'sessions'); continue; end
        ss = rd.(shanks{is}).sessions;
        for sIdx = 1:numel(ss)
            for ir = 1:2
                reg = regions{ir};
                if ~isfield(ss(sIdx), reg); continue; end
                rdata = ss(sIdx).(reg);
                if ~isfield(rdata, 'region_chan_depths_um'); continue; end
                dpt = double(rdata.region_chan_depths_um(:));
                if isempty(dpt); continue; end
                bounds.(reg)(1) = min(bounds.(reg)(1), min(dpt));
                bounds.(reg)(2) = max(bounds.(reg)(2), max(dpt));
            end
        end
    end
end
for ir = 1:2
    reg = regions{ir};
    if ~isfinite(bounds.(reg)(1)) || ~isfinite(bounds.(reg)(2))
        bounds.(reg) = [0, 6000];
    end
end
end


function T = build_channel_inclusion_long(SR_DATA, threshold_tag)
% One row per (rat, layer, region, session, window, included channel).
regions = {'PL', 'IL'};
windows = {'mono', 'poly'};
layers  = {'L2', 'L5'};
if strcmp(threshold_tag, '50pct'); chan_field = 'region_channels_50pct';
else;                              chan_field = 'region_channels_25pct';
end

rats = fieldnames(SR_DATA);
rats = rats(startsWith(rats, 'Rat'));
n_max  = numel(rats) * numel(layers) * 2 * numel(regions) * numel(windows);
small  = cell(n_max, 1);
n_used = 0;

for iR = 1:numel(rats)
    rd = SR_DATA.(rats{iR});
    if ~isfield(rd, 'L2_shank'); continue; end
    age_str      = string(rd.age);
    rat_str      = string(rd.rat_id);
    layer_shanks = [rd.L2_shank, rd.L5_shank];

    for li = 1:numel(layers)
        shFld = sprintf('Shank%d', layer_shanks(li));
        if ~isfield(rd, shFld); continue; end
        sess_struct = rd.(shFld).sessions;
        for sess = 1:numel(sess_struct)
            sess_name = string(sess_struct(sess).name);
            for ir = 1:numel(regions)
                if ~isfield(sess_struct(sess), regions{ir}); continue; end
                rdata = sess_struct(sess).(regions{ir});
                [chan_global, chan_depths] = get_region_chan_depths(rdata);
                if isempty(chan_global); continue; end

                for iw = 1:numel(windows)
                    if ~isfield(rdata, windows{iw}); continue; end
                    pk = rdata.(windows{iw});
                    if ~isstruct(pk) || ~isfield(pk, chan_field); continue; end
                    set_chans = pk.(chan_field)(:);
                    if isempty(set_chans); continue; end
                    in_set = ismember(chan_global, set_chans);
                    if ~any(in_set); continue; end

                    n_used = n_used + 1;
                    small{n_used} = make_subtable(rat_str, age_str, ...
                        string(layers{li}), string(regions{ir}), sess_name, ...
                        string(windows{iw}), chan_depths(in_set));
                end
            end
        end
    end
end

if n_used == 0; T = table(); return; end
T = vertcat(small{1:n_used});
end


function [chan_global, chan_depths] = get_region_chan_depths(rdata)
chan_global = []; chan_depths = [];
if isfield(rdata, 'region_chan_global_idx') && isfield(rdata, 'region_chan_depths_um')
    chan_global = rdata.region_chan_global_idx(:);
    chan_depths = double(rdata.region_chan_depths_um(:));
end
end


function T = make_subtable(rat, age, layer, reg, sess, win, depths)
n = numel(depths);
T = table(repmat(rat, n, 1), repmat(age, n, 1), repmat(layer, n, 1), ...
    repmat(reg, n, 1), repmat(sess, n, 1), repmat(win, n, 1), depths(:), ...
    'VariableNames', {'Rat','Age','Layer','Region','Session','Window','Depth_um'});
end
