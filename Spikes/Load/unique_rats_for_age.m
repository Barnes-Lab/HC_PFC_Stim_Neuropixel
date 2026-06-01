function rats = unique_rats_for_age(data, sess_types, neuron_type, age)
% Unique rat IDs that have any neurons in
% (age, neuron_type) across sess (NO INT FOR SOME SO THIS ACCOUNTS FOR
% THAT)


all_rats = {};
for si = 1:numel(sess_types)
    d = data.(sess_types{si}).(age).(neuron_type);
    if ~isempty(d.rat); all_rats = [all_rats; d.rat]; end
end
rats = unique(all_rats);
end
