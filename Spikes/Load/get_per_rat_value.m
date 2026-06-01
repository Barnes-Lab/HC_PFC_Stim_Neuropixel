function v = get_per_rat_value(data, sess_types, age, sess_filter, region_filter, ...
    win_filter, metric, thresh, rat_list)
% Per-rat  value for delta FR, recruitment, reliability (metrics) neural  (PYR ).
% per window (mono+poly) and  (PL+IL),
%

% SS 2026

% Window  mean of mono + poly per rat
if strcmp(win_filter, 'avg')
    v_m = get_per_rat_value(data, sess_types, age, sess_filter, region_filter, 'mono', metric, thresh, rat_list);
    v_p = get_per_rat_value(data, sess_types, age, sess_filter, region_filter, 'poly', metric, thresh, rat_list);
    v   = mean([v_m, v_p], 2, 'omitnan');
    return;
end

% ean of PL + IL per rat
if strcmp(region_filter, 'avg')
    v_PL = get_per_rat_value(data, sess_types, age, sess_filter, 'PL', win_filter, metric, thresh, rat_list);
    v_IL = get_per_rat_value(data, sess_types, age, sess_filter, 'IL', win_filter, metric, thresh, rat_list);
    v    = mean([v_PL, v_IL], 2, 'omitnan');
    return;
end

switch win_filter
    case 'mono'; wi = 1;
    case 'poly'; wi = 2;
    otherwise;   error('Unknown win_filter: %s', win_filter);
end
v = aggregate_neuron_weighted_per_rat(data, sess_types, 'PYR', age, ...
    sess_filter, region_filter, wi, metric, thresh, rat_list);
end
