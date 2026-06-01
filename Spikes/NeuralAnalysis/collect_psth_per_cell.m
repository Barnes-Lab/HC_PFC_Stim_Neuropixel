function res = collect_psth_per_cell(all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
    nt, age, region, sess_filter, time_bins, bsl_idx, PL_depth)
%PSTHs for (neuron_type, age, region)
% Pooled per-neuron dFR and per rat mean dFR. 
% Not using per rat score (amplified noise or biased toward neuron high count)
%
% %modified from Steinmentz github repo to include per rat and per cell
%SS 2026

neuron_psths = [];
neuron_rats  = {};

for u = 1:length(all_SP_filtered)
    if ~strcmp(neuron_classifications.type{u}, nt); continue; end
    if ~strcmp(all_SP_filtered(u).age, age);        continue; end

    if     isfield(all_SP_filtered(u), 'neuropixels_depth_uM'); depth = all_SP_filtered(u).neuropixels_depth_uM;
    elseif isfield(all_SP_filtered(u), 'depth_uM');             depth = all_SP_filtered(u).depth_uM;
    else;  continue;
    end
    nr = 'IL'; if depth < PL_depth; nr = 'PL'; end
    if ~strcmp(nr, region); continue; end

    rat = all_SP_filtered(u).rat;
    tEvents = [];
    for s = 1:length(sess_filter)
        T = ALL_STIM_TIMES(strcmp(string(ALL_STIM_TIMES.Rat), rat) & ...
                           strcmp(ALL_STIM_TIMES.Session, sess_filter{s}), :);
        if     ismember('Keep', T.Properties.VariableNames); T = T(~strcmp(T.Keep, 'N'), :);
        elseif ismember('keep', T.Properties.VariableNames); T = T(~strcmp(T.keep, 'N'), :);
        end
        if ~isempty(T); tEvents = [tEvents; T.Time_s]; end
    end
    if isempty(tEvents); continue; end

    spike_s = double(all_SP_filtered(u).t_uS) / 1e6;
    psth    = calculate_psth(spike_s, tEvents, time_bins);
    neuron_psths(end+1, :) = psth(:)';
    neuron_rats{end+1, 1}  = rat;
end

res.n_neurons = size(neuron_psths, 1);
if res.n_neurons == 0
    res.pooled_bsl = []; res.rat_bsl = []; res.n_rats = 0; return;
end

% Per-neuron baseline subtraction
neuron_bsl_mean = mean(neuron_psths(:, bsl_idx), 2);
neuron_bsl_subt = neuron_psths - neuron_bsl_mean;
res.pooled_bsl  = neuron_bsl_subt;

% Per-rat mean
unique_rats = unique(neuron_rats);
nR          = numel(unique_rats);
B           = size(neuron_psths, 2);
res.rat_bsl = nan(nR, B);
for r = 1:nR
    mask = strcmp(neuron_rats, unique_rats{r});
    res.rat_bsl(r, :) = mean(neuron_bsl_subt(mask, :), 1, 'omitnan');
end
res.n_rats = nR;
end
