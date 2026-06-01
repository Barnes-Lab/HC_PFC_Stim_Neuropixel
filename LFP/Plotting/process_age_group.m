function [rat_curves, avg_curve, sem_curve] = process_age_group( ...
    sub, MONO_DATA, stim_bins, bin_centers, amp_field)
%  Bin per-rat IL mono amplitudes into stim_bins and
% normalize each rat by subtracting its least-negative 
% Returns per-rat curves plus across-rat mean and SEM.
%
%
% SS 2026

n_rats     = height(sub);
n_bins     = numel(bin_centers);
rat_curves = nan(n_rats, n_bins);

for i = 1:n_rats
    ratFld = sprintf('Rat%s',   sub.Rat{i});
    shFld  = sprintf('Shank%d', sub.Shank(i));
    sess   = sub.Session(i);

    if ~isfield(MONO_DATA, ratFld) || ~isfield(MONO_DATA.(ratFld), shFld); continue; end
    sessions = MONO_DATA.(ratFld).(shFld).sessions;
    if sess > numel(sessions); continue; end
    sd = sessions(sess);
    if ~isfield(sd, 'IL') || ~isstruct(sd.IL.mono); continue; end
    if ~isfield(sd.IL.mono, amp_field); continue; end

    currents = sd.stim_currents(sd.valid_stims(:));
    amps     = sd.IL.mono.(amp_field)(:);
    if numel(currents) ~= numel(amps); continue; end

    keep = ~isnan(currents) & ~isnan(amps);
    currents = currents(keep); amps = amps(keep);
    if isempty(currents); continue; end

    binned = nan(1, n_bins);
    for b = 1:n_bins
        in = currents >= stim_bins(b) & currents < stim_bins(b+1);
        if any(in); binned(b) = mean(amps(in)); end
    end

    base_amp = max(binned);
    if isnan(base_amp); continue; end
    rat_curves(i, :) = binned - base_amp;
end

avg_curve = mean(rat_curves, 1, 'omitnan');
sem_curve = std(rat_curves, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(rat_curves), 1));
end
