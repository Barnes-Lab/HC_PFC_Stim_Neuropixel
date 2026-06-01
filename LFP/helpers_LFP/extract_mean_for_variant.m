function [mean_amp, mean_lat] = extract_mean_for_variant( ...
    win_data, amp_field, lat_field)
% Mean across the  5 analysis stims 
mean_amp = NaN;
mean_lat = NaN;
if ~isfield(win_data, amp_field) || ~isfield(win_data, lat_field)
    return;
end
amps = win_data.(amp_field);
lats = win_data.(lat_field);
if any(~isnan(amps)); mean_amp = nanmean(amps); end
if any(~isnan(lats)); mean_lat = nanmean(lats); end
end
