function neuron_classifications = classify_neurons_filtered(all_SP_filtered, baseDir, apply_artifact_rejection)
%   Classify filtered units into INT / PYR /UNCLASSIFIED 
% based on waveform trough-to-peak and half-width thresholds,
% Saves classifications.mat and
% a classification summary plot 
% DIR/neuron_classification/.
%
% Thresholds (From Bartho 2004):
%   INT  trough-to-peak-[0.16, 0.50] ms && half-width < 0.5 ms
%   PYR  trough-to-peak [0.50, 1.13] ms
%
% OUTPUT
%   neuron_classifications  struct with per-unit type, waveform, confidence,
%                           trough_to_peak, half_width, ac_peak_latency,
%                           ac_decay_log_slope, artifact_free
%
% SS 2026

if nargin < 3; apply_artifact_rejection = false; end

fs = 30000;
thresholds = struct( ...
    'trough_to_peak_int_lower', 0.16, ...
    'trough_to_peak_int_upper', 0.50, ...
    'trough_to_peak_pyr_lower', 0.50, ...
    'trough_to_peak_pyr_upper', 1.13, ...
    'half_width_narrow',        0.5);

n_units = length(all_SP_filtered);
neuron_classifications = struct();
neuron_classifications.type               = cell(n_units, 1);
neuron_classifications.waveform           = cell(n_units, 1);
neuron_classifications.confidence         = zeros(n_units, 1);
neuron_classifications.trough_to_peak     = zeros(n_units, 1);
neuron_classifications.half_width         = zeros(n_units, 1);
neuron_classifications.ac_peak_latency    = zeros(n_units, 1);
neuron_classifications.ac_decay_log_slope = zeros(n_units, 1);
neuron_classifications.artifact_free      = true(n_units, 1);

figDir = fullfile(baseDir, 'neuron_classification');
if ~exist(figDir, 'dir'); mkdir(figDir); end
fprintf('Classifying %d neurons...\n', n_units);

for u = 1:n_units
    SP = all_SP_filtered(u);
    if ~isfield(SP, 'WV') || ~isfield(SP.WV, 'mWV')
        neuron_classifications.type{u}                = 'UNCLASSIFIED';
        neuron_classifications.waveform{u}            = 'UNCLASSIFIED';
        neuron_classifications.confidence(u)          = 0;
        neuron_classifications.ac_peak_latency(u)     = NaN;
        neuron_classifications.ac_decay_log_slope(u)  = NaN;
        continue;
    end

    wf_params = get_waveform_spec(SP.WV.mWV(:, 1)', fs, false);
    if isfield(SP, 't_uS')
        ac_params = get_autocorr_features(SP.t_uS, 1, 500, false);
    else
        ac_params = struct('classification', 'UNCLASSIFIED', ...
            'peak_latency_ms', NaN, 'decay_log_slope', NaN);
    end

    neuron_classifications.trough_to_peak(u)     = wf_params.trough_to_peak_ms;
    neuron_classifications.half_width(u)         = wf_params.half_width_ms;
    neuron_classifications.ac_peak_latency(u)    = ac_params.peak_latency_ms;
    neuron_classifications.ac_decay_log_slope(u) = ac_params.decay_log_slope;

    tp = wf_params.trough_to_peak_ms;
    hw = wf_params.half_width_ms;
    if ~isnan(tp) && ~isnan(hw)
        if     tp >= thresholds.trough_to_peak_int_lower && ...
               tp <= thresholds.trough_to_peak_int_upper && ...
               hw <  thresholds.half_width_narrow
            neuron_classifications.type{u}       = 'INT';
            neuron_classifications.confidence(u) = 0.9;
        elseif tp >= thresholds.trough_to_peak_pyr_lower && ...
               tp <= thresholds.trough_to_peak_pyr_upper
            neuron_classifications.type{u}       = 'PYR';
            neuron_classifications.confidence(u) = 0.9;
        else
            neuron_classifications.type{u}       = 'UNCLASSIFIED';
            neuron_classifications.confidence(u) = 0.3;
        end

        switch ac_params.classification
            case 'FS'; neuron_classifications.waveform{u} = 'FS';
            case 'BT'; neuron_classifications.waveform{u} = 'BS';   % BT -> BS relabel
            case 'RF'; neuron_classifications.waveform{u} = 'RS';
            otherwise; neuron_classifications.waveform{u} = 'UNCLASSIFIED';
        end
    else
        neuron_classifications.type{u}       = 'UNCLASSIFIED';
        neuron_classifications.waveform{u}   = 'UNCLASSIFIED';
        neuron_classifications.confidence(u) = 0;
    end
end

if apply_artifact_rejection
    fprintf('\nApplying artifact rejection...\n');
    neuron_classifications = apply_artifact_filter(neuron_classifications, all_SP_filtered, figDir);
end

save(fullfile(figDir, 'neuron_classifications.mat'), 'neuron_classifications');

%% Summary
INT_total     = sum(strcmp(neuron_classifications.type, 'INT'));
PYR_total     = sum(strcmp(neuron_classifications.type, 'PYR'));
UNCLASS_total = sum(strcmp(neuron_classifications.type, 'UNCLASSIFIED'));
fprintf('\nClassification Summary:\n');
fprintf('INT: %d neurons\nPYR: %d neurons\nUnclassified: %d neurons\n', ...
    INT_total, PYR_total, UNCLASS_total);

if apply_artifact_rejection
    INT_clean = sum(strcmp(neuron_classifications.type, 'INT') & neuron_classifications.artifact_free);
    PYR_clean = sum(strcmp(neuron_classifications.type, 'PYR') & neuron_classifications.artifact_free);
    fprintf('\nAfter artifact rejection:\n');
    fprintf('INT: %d/%d (%.1f%% retained)\n', INT_clean, INT_total, 100*INT_clean/INT_total);
    fprintf('PYR: %d/%d (%.1f%% retained)\n', PYR_clean, PYR_total, 100*PYR_clean/PYR_total);
end

create_classification_summary_plot(neuron_classifications, all_SP_filtered, figDir, apply_artifact_rejection);
end
