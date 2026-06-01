function plot_within_rat_panel(ax, vY_a, vY_b, vO_a, vO_b, ageCols, ...
    xtick_labels, ylab, ttl, varargin)
% PLOT_WITHIN_RAT_PANEL  One within-rat paired panel: 2 conditions on x-axis,
% individual rats as light lines (background, both ages), group mean line
% (LineWidth 5) + SEM errorbars on top so neither age's rats overlay the
% other group's mean.
%
% Example:
%   plot_within_rat_panel(gca, vY_a, vY_b, vO_a, vO_b, ageCols, ...
%       {'Mono','Poly'}, '\Delta FR (Hz)', 'PL', 'PRM_PUBIFY_SIZE', 16);
%
% INPUT
%   ax            axes handle
%   vY_a, vY_b    [n_rats_Y x 1]  paired young values
%   vO_a, vO_b    [n_rats_O x 1]  paired old values
%   ageCols       struct .Y .O
%   xtick_labels  {1x2}
%   ylab, ttl     char
%   varargin      'PRM_PUBIFY_SIZE' (default 14)
%
% SS 2026

PRM_PUBIFY_SIZE = 14;
Extract_varargin;

axes(ax); hold(ax, 'on');

[Y_a, Y_b, nY] = clean_pair(vY_a, vY_b);
[O_a, O_b, nO] = clean_pair(vO_a, vO_b);

draw_rats(ax, Y_a, Y_b, ageCols.Y);
draw_rats(ax, O_a, O_b, ageCols.O);
draw_mean(ax, Y_a, Y_b, ageCols.Y, nY);
draw_mean(ax, O_a, O_b, ageCols.O, nO);

xlim(ax, [0.5, 2.5]);
set(ax, 'XTick', [1, 2], 'XTickLabel', xtick_labels);
ylabel(ax, ylab);
title(ax, ttl);
box(ax, 'off');
if exist('pubify_figure_axis_robust', 'file')
    pubify_figure_axis_robust(PRM_PUBIFY_SIZE, PRM_PUBIFY_SIZE);
end
end


function [a, b, n] = clean_pair(va, vb)
ok = ~isnan(va) & ~isnan(vb);
a  = va(ok); b = vb(ok); n = numel(a);
end


function draw_rats(ax, vA, vB, col)
if isempty(vA); return; end
light = col * 0.3 + [1 1 1] * 0.7;
for r = 1:numel(vA)
    plot(ax, [1, 2], [vA(r), vB(r)], 'o-', 'Color', light, 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', light, 'HandleVisibility', 'off');
end
end


function draw_mean(ax, vA, vB, col, n)
if n == 0; return; end
mA = mean(vA); sA = std(vA) / sqrt(n);
mB = mean(vB); sB = std(vB) / sqrt(n);
plot(ax, [1, 2], [mA, mB], '-', 'Color', col, 'LineWidth', 5);
errorbar(ax, 1, mA, sA, 'Color', col, 'LineWidth', 2.5, 'CapSize', 10, 'HandleVisibility', 'off');
errorbar(ax, 2, mB, sB, 'Color', col, 'LineWidth', 2.5, 'CapSize', 10, 'HandleVisibility', 'off');
end
