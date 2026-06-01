function plot_bar_with_sem(ax, v, xc, col, varargin)
% Mean bar + SEM errorbar at 1 x.+
% SS 2026

PRM_WIDTH       = 0.65;
PRM_FACE_ALPHA  = 0.5;
PRM_EDGE        = 'col';
PRM_LINE_WIDTH  = 1.5;
PRM_ERR_LW      = 2;
PRM_CAP_SIZE    = 8;
Extract_varargin;

if isempty(v); return; end

if     ischar(PRM_EDGE) && strcmp(PRM_EDGE, 'col'); edge_col = col;
elseif ischar(PRM_EDGE) && strcmp(PRM_EDGE, 'k');   edge_col = 'k';
else;                                                edge_col = PRM_EDGE;
end

m  = mean(v, 'omitnan');
se = std(v, 'omitnan') / sqrt(sum(~isnan(v)));

bar(ax, xc, m, PRM_WIDTH, 'FaceColor', col, 'FaceAlpha', PRM_FACE_ALPHA, ...
    'EdgeColor', edge_col, 'LineWidth', PRM_LINE_WIDTH);
errorbar(ax, xc, m, se, 'Color', edge_col, 'LineStyle', 'none', ...
    'LineWidth', PRM_ERR_LW, 'CapSize', PRM_CAP_SIZE, 'HandleVisibility', 'off');
end
