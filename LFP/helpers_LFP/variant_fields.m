function [amp_field, lat_field] = variant_fields(variant)
% Map variant -> SR_DATA peak-struct field names.
switch lower(variant)
    case 'lowest'
        amp_field = 'best_channel_amplitude';
        lat_field = 'best_channel_latency';
    case '3chan'
        amp_field = 'mean_amplitude_3chan';
        lat_field = 'mean_latency_3chan';
    case '50pct'
        amp_field = 'mean_amplitude_50pct';
        lat_field = 'mean_latency_50pct';
    case '25pct'
        amp_field = 'mean_amplitude_25pct';
        lat_field = 'mean_latency_25pct';
    otherwise
        error(['Unknown amplitude_variant ''%s''. ' ...
               'Use lowest | 3chan | 50pct | 25pct'], variant);
end
end
