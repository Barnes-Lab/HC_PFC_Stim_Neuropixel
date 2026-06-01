function plot_scatter_with_jitter(ax, v, xc, col, varargin)
% CHECK FOR PARAMS BELOW
PRM_MARKER_SIZE = 50;
PRM_JITTER      = 0.12;
PRM_LIGHTEN     = true;
PRM_EDGE        = 'auto';
PRM_FACE_ALPHA  = 1;
Extract_varargin;

if isempty(v); return; end

if PRM_LIGHTEN; face_col = col*0.3 + [1 1 1]*0.7;
else;           face_col = col;
end
if     ischar(PRM_EDGE) && strcmp(PRM_EDGE, 'auto'); edge_col = col*0.5 + [1 1 1]*0.5;
elseif ischar(PRM_EDGE) && strcmp(PRM_EDGE, 'none'); edge_col = 'none';
else;                                                edge_col = PRM_EDGE;
end

jit = (rand(numel(v), 1) - 0.5) * PRM_JITTER;
scatter(ax, ones(numel(v), 1)*xc + jit, v(:), PRM_MARKER_SIZE, face_col, 'filled', ...
    'MarkerEdgeColor', edge_col, 'MarkerFaceAlpha', PRM_FACE_ALPHA, ...
    'HandleVisibility', 'off');
end
