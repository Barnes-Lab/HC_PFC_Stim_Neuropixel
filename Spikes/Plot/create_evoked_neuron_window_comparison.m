function create_evoked_neuron_window_comparison(data, sess_types, neuron_types, ageCols, figDir)
% per neuron and per rat comparisons per window and reigon and combinations
% of all
%
% SS 2026

EVOKED_THR = 0;
regions    = {'PL', 'IL'};
ages       = {'Y', 'O'};

cat_colors_full = [0.30 0.50 0.80;   % Mono
                   0.80 0.50 0.30;   % Poly
                   0.20 0.70 0.20;   % Both
                   0.70 0.70 0.70];  % Neither
cat_colors_resp = cat_colors_full(1:3, :);
cats_full = {'Mono', 'Poly', 'Both', 'Neither'};
cats_resp = {'Mono', 'Poly', 'Both'};

result = struct();
for ri = 1:2
    for ai = 1:2
        for ti = 1:numel(neuron_types)
            result.(regions{ri}).(ages{ai}).(neuron_types{ti}) = ...
                collect_cell(data, sess_types, regions{ri}, ages{ai}, neuron_types{ti}, EVOKED_THR);
        end
    end
end

plot_pies('pooled', 'all',        result, regions, ages, neuron_types, cats_full, cat_colors_full, figDir);
plot_pies('pooled', 'responders', result, regions, ages, neuron_types, cats_resp, cat_colors_resp, figDir);
plot_pies('perRat', 'all',        result, regions, ages, neuron_types, cats_full, cat_colors_full, figDir);
plot_pies('perRat', 'responders', result, regions, ages, neuron_types, cats_resp, cat_colors_resp, figDir);

plot_evoked_FR(result, regions, ages, neuron_types, ageCols, figDir);
end


function res = collect_cell(data, sess_types, region, age, type, EVOKED_THR)
w1 = []; w2 = []; rats = {};
for si = 1:numel(sess_types)
    d = data.(sess_types{si}).(age).(type);
    if isempty(d.rat); continue; end
    rmask = strcmp(d.region, region);
    if ~any(rmask); continue; end
    w1   = [w1;   d.win1(rmask)];
    w2   = [w2;   d.win2(rmask)];
    rats = [rats; d.rat(rmask)];
end

res = struct('n_total', 0, 'n_mono', 0, 'n_poly', 0, 'n_both', 0, 'n_neither', 0, ...
    'mono_fr', [], 'poly_fr', [], ...
    'rat_pct_mono', [], 'rat_pct_poly', [], 'rat_pct_both', [], 'rat_pct_neither', [], ...
    'rat_n', []);

res.n_total = numel(w1);
if res.n_total == 0; return; end

im = w1 > EVOKED_THR; ip = w2 > EVOKED_THR;
res.n_mono    = sum(im & ~ip);
res.n_poly    = sum(~im & ip);
res.n_both    = sum(im & ip);
res.n_neither = sum(~im & ~ip);
res.mono_fr   = w1(im);
res.poly_fr   = w2(ip);

unique_rats = unique(rats);
nR = numel(unique_rats);
res.rat_pct_mono    = zeros(nR, 1);
res.rat_pct_poly    = zeros(nR, 1);
res.rat_pct_both    = zeros(nR, 1);
res.rat_pct_neither = zeros(nR, 1);
res.rat_n           = zeros(nR, 1);
for r = 1:nR
    rm   = strcmp(rats, unique_rats{r});
    nr   = sum(rm);
    im_r = im(rm); ip_r = ip(rm);
    res.rat_pct_mono(r)    = 100 * sum(im_r & ~ip_r) / nr;
    res.rat_pct_poly(r)    = 100 * sum(~im_r & ip_r) / nr;
    res.rat_pct_both(r)    = 100 * sum(im_r & ip_r) / nr;
    res.rat_pct_neither(r) = 100 * sum(~im_r & ~ip_r) / nr;
    res.rat_n(r)           = nr;
end
end


function plot_pies(scope, mode, result, regions, ages, neuron_types, cats, cat_colors, figDir)
is_resp = strcmp(mode, 'responders');
nT      = numel(neuron_types);
nRows   = 2 * nT;

figure('Position', [100, 100, 560, 180*nRows]);
tiledlayout(nRows, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

tile_i = 0;
for ri = 1:numel(regions)
    for ti = 1:nT
        for ai = 1:2
            tile_i = tile_i + 1;
            ax  = nexttile;
            res = result.(regions{ri}).(ages{ai}).(neuron_types{ti});
            [vals, n_disp] = compute_slice_values(res, scope, is_resp);
            if isempty(vals) || sum(vals) == 0
                axis(ax, 'off');
                title(ax, sprintf('%s-%s %s', regions{ri}, neuron_types{ti}, ages{ai}));
                continue;
            end
            h = pie(ax, vals);
            apply_pie_colors(h, vals, cat_colors);
            title(ax, sprintf('%s-%s %s (n=%d)', regions{ri}, neuron_types{ti}, ages{ai}, n_disp));
            set(ax, 'FontSize', 9);
            if tile_i == 1
                lg = legend(cats, 'Location', 'eastoutside', 'FontSize', 9);
                lg.Layout.Tile = 'east';
            end
        end
    end
end

bias_str = ''; if strcmp(scope, 'pooled'); bias_str = ' (biased by N)'; end
sgtitle(sprintf('Evoked Neurons - %s, %s%s', mode, scope, bias_str), 'FontSize', 12);
saveas(gcf, fullfile(figDir, sprintf('evoked_neuron_proportions_%s_%s.png', scope, mode)));
close(gcf);
end


function [vals, n_disp] = compute_slice_values(res, scope, is_resp)
if res.n_total == 0; vals = []; n_disp = 0; return; end

switch scope
    case 'pooled'
        if is_resp
            n_use = res.n_mono + res.n_poly + res.n_both;
            if n_use == 0; vals = []; n_disp = 0; return; end
            vals = [res.n_mono, res.n_poly, res.n_both];
            n_disp = n_use;
        else
            vals = [res.n_mono, res.n_poly, res.n_both, res.n_neither];
            n_disp = res.n_total;
        end
    case 'perRat'
        nR = numel(res.rat_pct_mono);
        if nR == 0; vals = []; n_disp = 0; return; end
        if is_resp
            tot = res.rat_pct_mono + res.rat_pct_poly + res.rat_pct_both;
            ok  = tot > 0;
            if ~any(ok); vals = []; n_disp = 0; return; end
            pm = res.rat_pct_mono(ok) ./ tot(ok);
            pp = res.rat_pct_poly(ok) ./ tot(ok);
            pb = res.rat_pct_both(ok) ./ tot(ok);
            vals = [mean(pm), mean(pp), mean(pb)];
            n_disp = sum(ok);
        else
            vals = [mean(res.rat_pct_mono), mean(res.rat_pct_poly), ...
                    mean(res.rat_pct_both), mean(res.rat_pct_neither)];
            n_disp = nR;
        end
end
end


function apply_pie_colors(h, vals, cat_colors)
patch_h = h(arrayfun(@(x) isa(x, 'matlab.graphics.primitive.Patch'), h));
nz = find(vals > 0);
for k = 1:numel(patch_h)
    if k <= numel(nz)
        patch_h(k).FaceColor = cat_colors(nz(k), :);
        patch_h(k).EdgeColor = 'k';
        patch_h(k).LineWidth = 0.8;
    end
end
end


function plot_evoked_FR(result, regions, ages, neuron_types, ageCols, figDir)
% Mono = lighter alpha (0.40); Poly = stronger alpha (0.85).
% Layering: scatter all neurons first (background), then bars+errorbars on top.
nT = numel(neuron_types);

all_fr = [];
for ri = 1:2
    for ai = 1:2
        for ti = 1:nT
            r = result.(regions{ri}).(ages{ai}).(neuron_types{ti});
            all_fr = [all_fr; r.mono_fr(:); r.poly_fr(:)];
        end
    end
end
all_fr = all_fr(isfinite(all_fr));
if exist('get_unified_ylim', 'file');  fr_yl = get_unified_ylim({all_fr}, 0.20);
elseif isempty(all_fr);                fr_yl = [-1, 1];
else;                                   fr_yl = [min(all_fr) - 0.1*range(all_fr), max(all_fr) + 0.1*range(all_fr)];
end

x_map = struct('Y', [1, 2], 'O', [4, 5]);

figure('Position', [100, 100, 300*nT, 500]);
tiledlayout(2, nT, 'TileSpacing', 'compact', 'Padding', 'compact');
for ri = 1:2
    for ti = 1:nT
        ax = nexttile; hold(ax, 'on');

        % Pass 1: small no-edge scatter for all neurons
        for ai = 1:2
            res = result.(regions{ri}).(ages{ai}).(neuron_types{ti});
            xp  = x_map.(ages{ai});
            plot_scatter_with_jitter(ax, res.mono_fr, xp(1), ageCols.(ages{ai}), ...
                'PRM_MARKER_SIZE', 18, 'PRM_JITTER', 0.1, 'PRM_EDGE', 'none', 'PRM_FACE_ALPHA', 0.5);
            plot_scatter_with_jitter(ax, res.poly_fr, xp(2), ageCols.(ages{ai}), ...
                'PRM_MARKER_SIZE', 18, 'PRM_JITTER', 0.1, 'PRM_EDGE', 'none', 'PRM_FACE_ALPHA', 0.5);
        end

        % Pass 2: bars with k edge, variable face alpha
        for ai = 1:2
            res = result.(regions{ri}).(ages{ai}).(neuron_types{ti});
            xp  = x_map.(ages{ai});
            plot_bar_with_sem(ax, res.mono_fr, xp(1), ageCols.(ages{ai}), ...
                'PRM_WIDTH', 0.7, 'PRM_FACE_ALPHA', 0.40, 'PRM_EDGE', 'k', ...
                'PRM_LINE_WIDTH', 1.2, 'PRM_ERR_LW', 1.2);
            plot_bar_with_sem(ax, res.poly_fr, xp(2), ageCols.(ages{ai}), ...
                'PRM_WIDTH', 0.7, 'PRM_FACE_ALPHA', 0.85, 'PRM_EDGE', 'k', ...
                'PRM_LINE_WIDTH', 1.2, 'PRM_ERR_LW', 1.2);
        end

        plot(ax, [3, 3], fr_yl, 'k--', 'LineWidth', 1);
        xlim(ax, [0.4, 5.6]); ylim(ax, fr_yl);
        set(ax, 'XTick', [1.5, 4.5], 'XTickLabel', {'Y', 'O'});
        ylabel(ax, '\Delta FR (Hz)');
        title(ax, sprintf('%s - %s', regions{ri}, neuron_types{ti}));
        box(ax, 'off');
        if exist('pubify_figure_axis_robust', 'file'); pubify_figure_axis_robust(12, 12); end
    end
end
sgtitle('\Delta FR of Evoked Neurons (Mono = light, Poly = solid)');
saveas(gcf, fullfile(figDir, 'evoked_neuron_FR_mono_vs_poly.png'));
close(gcf);
end
