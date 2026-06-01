function [u_start, u_end] = union_crossings(out_a, out_b, fld_start, fld_end)
% Union: earliest start across valid layers, latest end across valid layers.
% Skip a layer if its .valid == false. If both invalid, return NaN/NaN.
starts = [];
ends   = [];
if out_a.valid; starts = [starts, out_a.(fld_start)]; ends = [ends, out_a.(fld_end)]; end
if out_b.valid; starts = [starts, out_b.(fld_start)]; ends = [ends, out_b.(fld_end)]; end

if isempty(starts)
    u_start = NaN; u_end = NaN;
else
    u_start = min(starts);
    u_end   = max(ends);
end
end
