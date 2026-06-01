function v = aggregate_neuron_weighted_per_rat(data, sess_types, neuron_type, age, ...
    sess_filter, region_filter, win_idx, metric, thresh, rat_list)
%
%Pool all neurons per ses/region per rat and aggreate metric 
% (deltaFR |Recruitment | Reliabttility)
%  NaN if rat has no matching neurons.
%

if win_idx == 1; w_f = 'win1'; r_f = 'win1_reliability';
else;            w_f = 'win2'; r_f = 'win2_reliability';
end

if strcmp(sess_filter, 'all'); sl = sess_types;
else;                          sl = {sess_filter};
end

n = numel(rat_list);
v = nan(n, 1);

for r = 1:n
    rat = rat_list{r};
    n_total = count_matching(data, sl, age, neuron_type, rat, region_filter);
    if n_total == 0; continue; end

    pw = nan(n_total, 1);
    pr = nan(n_total, 1);
    ofs = 0;
    for si = 1:numel(sl)
        d = data.(sl{si}).(age).(neuron_type);
        if isempty(d.rat); continue; end
        m = strcmp(d.rat, rat);
        if ~strcmp(region_filter, 'all'); m = m & strcmp(d.region, region_filter); end
        n_si = sum(m);
        if n_si == 0; continue; end
        pw(ofs+1 : ofs+n_si) = d.(w_f)(m);
        pr(ofs+1 : ofs+n_si) = d.(r_f)(m);
        ofs = ofs + n_si;
    end

    switch metric
        case 'delta_FR';              v(r) = mean(pw, 'omitnan');
        case 'recruitment';           v(r) = 100 * sum(pw > thresh) / n_total;
        case 'reliability_per_trial'; v(r) = mean(pr, 'omitnan');
        otherwise;                    error('Unknown metric: %s', metric);
    end
end
end


function n_total = count_matching(data, sl, age, neuron_type, rat, region_filter)
n_total = 0;
for si = 1:numel(sl)
    d = data.(sl{si}).(age).(neuron_type);
    if isempty(d.rat); continue; end
    m = strcmp(d.rat, rat);
    if ~strcmp(region_filter, 'all'); m = m & strcmp(d.region, region_filter); end
    n_total = n_total + sum(m);
end
end
