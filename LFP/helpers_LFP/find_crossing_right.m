function t_cross = find_crossing_right(trace, t_ms, peak_idx, thresh)
% Walk right from peak_idx until trace >= thresh.
n = numel(trace);
t_cross = t_ms(end);
for k = peak_idx:(n-1)
    if trace(k) < thresh && trace(k+1) >= thresh
        t_cross = interp_cross(t_ms(k), t_ms(k+1), trace(k), trace(k+1), thresh);
        return;
    end
end
end
