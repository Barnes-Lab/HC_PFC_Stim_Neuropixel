function params = get_waveform_spec(mean_wf, fs, plot_it, varargin)
% Extract waveform parameters (trough-to-peak,
% half-width, AHP) with spline upsampling. Baseline-subtracted using
% samples 15-25.
%
% INPUT
%   mean_wf
%   fs  Hzdefault 30000)
%   plot_it  default false (takes too long otherwise)
%
%
% SS 2026

if nargin < 2; fs = 30000; end
if nargin < 3; plot_it = false; end
n_up_pts = 1000;

rat     = '10947';
unit_id = 1;
Extract_varargin;

fig_dir = fullfile('C:\DATA\NP2', 'waveform_analysis', rat);
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

% Baseline-subtract using samples 15-25 % COWEN 
mean_wf = mean_wf - mean(mean_wf(15:25));

% Time
tmp           = (1:length(mean_wf)) - 1;
wv_x_ms       = 1000 * (tmp / fs);
wv_x_ms_up    = linspace(wv_x_ms(1), wv_x_ms(end), n_up_pts);
up_ms_per_samp = (wv_x_ms_up(end) - wv_x_ms_up(1)) / n_up_pts;

% Upsample and normalize to trough (SAME AS COWEN)
wave_up   = interp1(wv_x_ms, mean_wf, wv_x_ms_up, 'spline');
wave_norm = wave_up / abs(min(wave_up));

% Trough AHPC
[trough_val, trough_idx] = min(wave_norm);
wave_after_trough = wave_norm;
wave_after_trough(1:trough_idx) = -Inf;
[peak_val, peak_idx] = max(wave_after_trough);

params.peak_threshold_used = 0;
params.trough_to_peak_ms   = (peak_idx - trough_idx) * up_ms_per_samp;

% Half-width 
half_depth = trough_val / 2;
crossings  = find(wave_norm <= half_depth);
if ~isempty(crossings)
    params.half_width_ms       = (crossings(end) - crossings(1)) * up_ms_per_samp;
    params.half_width_start_idx = crossings(1);
    params.half_width_end_idx   = crossings(end);
else
    params.half_width_ms        = NaN;
    params.half_width_start_idx = NaN;
    params.half_width_end_idx   = NaN;
end

% AHP
wave_after_peak = wave_norm;
wave_after_peak(1:peak_idx) = 0;
params.AHP = sum(wave_after_peak(wave_after_peak > 0));

params.trough_val     = trough_val;
params.trough_idx     = trough_idx;
params.peak_val       = peak_val;
params.peak_idx       = peak_idx;
params.wave_norm      = wave_norm;
params.wv_x_ms_up     = wv_x_ms_up;
params.up_ms_per_samp = up_ms_per_samp;

if plot_it
    wv_x_ms_aligned = wv_x_ms_up - wv_x_ms_up(trough_idx);
    figure('Position', [50 50 1100 800], 'Visible', 'off');
    plot(wv_x_ms_aligned, wave_norm, 'k', 'LineWidth', 2); hold on;
    plot(0, trough_val, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    plot(wv_x_ms_aligned(peak_idx), peak_val, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    plot([0, wv_x_ms_aligned(peak_idx)], [trough_val, trough_val], 'r-', 'LineWidth', 2);
    plot([wv_x_ms_aligned(peak_idx), wv_x_ms_aligned(peak_idx)], [-1, peak_val], 'k--', 'LineWidth', 1);
    if ~isnan(params.half_width_ms)
        plot([wv_x_ms_aligned(params.half_width_start_idx), wv_x_ms_aligned(params.half_width_end_idx)], ...
            [half_depth, half_depth], 'g-', 'LineWidth', 2);
    end
    xlabel('Time from Trough (ms)');
    ylabel('Normalized Amplitude');
    title(sprintf('Trough-to-Peak: %.3fms, Half-Width: %.3fms', ...
        params.trough_to_peak_ms, params.half_width_ms));
    xlim([-1 1]);
    pubify_figure_axis_robust(15, 15);
    box off;
    saveas(gcf, fullfile(fig_dir, sprintf('%03d_waveform_params.png', unit_id)));
    close;
end
end
