function tc = interp_cross(t1, t2, v1, v2, thresh)
% Linear interpolation: when does the line through (t1,v1)-(t2,v2) hit thresh?
if v2 == v1; tc = t1; return; end
tc = t1 + (thresh - v1) * (t2 - t1) / (v2 - v1);
end
