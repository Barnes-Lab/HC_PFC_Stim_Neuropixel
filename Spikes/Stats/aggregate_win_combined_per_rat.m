function v = aggregate_win_combined_per_rat(data, sess_types, age, ...
    sess_filter, region_filter, responders_only, recruit_thresh, rat_list)
% AGGREGATE_WIN_COMBINED_PER_RAT  Per-rat mean of win_combined (PYR only).
% Optionally restricts to responders (win1 > thresh OR win2 > thresh).
%
% SS 2026

if strcmp(sess_filter, 'all'); sl = sess_types;
else;                          sl = {sess_filter};
end

n = numel(rat_list);
v = nan(n, 1);

for r = 1:n
    rat = rat_list{r};
    n_total = 0;
    for si = 1:numel(sl)
        d = data.(sl{si}).(age).PYR;
        if isempty(d.rat); continue; end
        m = strcmp(d.rat, rat);
        if ~strcmp(region_filter, 'all'); m = m & strcmp(d.region, region_filter); end
        n_total = n_total + sum(m);
    end
    if n_total == 0; continue; end

    wc = nan(n_total, 1); w1 = nan(n_total, 1); w2 = nan(n_total, 1);
    ofs = 0;
    for si = 1:numel(sl)
        d = data.(sl{si}).(age).PYR;
        if isempty(d.rat); continue; end
        m = strcmp(d.rat, rat);
        if ~strcmp(region_filter, 'all'); m = m & strcmp(d.region, region_filter); end
        k = sum(m);
        if k == 0; continue; end
        wc(ofs+1 : ofs+k) = d.win_combined(m);
        w1(ofs+1 : ofs+k) = d.win1(m);
        w2(ofs+1 : ofs+k) = d.win2(m);
        ofs = ofs + k;
    end

    if responders_only
        wc = wc((w1 > recruit_thresh) | (w2 > recruit_thresh));
    end
    if isempty(wc) || all(isnan(wc)); continue; end
    v(r) = mean(wc, 'omitnan');
end
end
