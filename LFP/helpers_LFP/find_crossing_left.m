function t_cross = find_crossing_left(trace, t_ms, peak_idx, thresh)

t_cross = t_ms(1);
for k = peak_idx:-1:2
    if trace(k-1) >= thresh && trace(k) < thresh
        % Linear interpolation between samples k-1 and k
        t_cross = interp_cross(t_ms(k-1), t_ms(k), trace(k-1), trace(k), thresh);
        return;
    end
end
end
