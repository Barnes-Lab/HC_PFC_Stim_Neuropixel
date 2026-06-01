function ac_params = get_autocorr_features(spike_times_us, bin_size_ms, max_lag_ms, plot_it)
% Autocorrelation-based classification  Computes ISI histogram, peak latency, log-slope
% of exponential decay after peak, 
%
% DID NOT USE - FR VERY DIFFERENT UNDER ANESTHESIA
%
% SS 2026

if nargin < 2; bin_size_ms = 1;   end
if nargin < 3; max_lag_ms  = 500; end
if nargin < 4; plot_it     = false; end

%% Autocorrelogram via ISI histogram
spike_times_s = spike_times_us / 1e6;
max_lag_s     = max_lag_ms / 1000;
bin_size_s    = bin_size_ms / 1000;
edges         = -max_lag_s : bin_size_s : max_lag_s;
lags_s        = edges(1:end-1) + bin_size_s/2;

isi = diff(spike_times_s);
ac  = histcounts(isi, edges);
ac  = ac / max(ac);

%% Peak latency (post-zero-lag)
zero_idx     = find(lags_s >= 0, 1, 'first');
search_range = zero_idx:length(ac);
[~, peak_idx_rel] = max(ac(search_range));
peak_idx        = search_range(1) + peak_idx_rel - 1;
peak_latency_ms = lags_s(peak_idx) * 1000;

%% Exponential decay fit on log(ac) after peak
fit_range = peak_idx : min(peak_idx + 100, length(ac));
if length(fit_range) > 10
    x     = lags_s(fit_range) * 1000;
    y     = ac(fit_range);
    valid = y > 0;
    if sum(valid) > 5
        p           = polyfit(x(valid), log(y(valid)), 1);
        decay_slope = p(1);
        decay_log_slope = log10(abs(decay_slope));
    else
        decay_log_slope = NaN;
    end
else
    decay_log_slope = NaN;
end

%% Classification (Xing et al. 2020 criteria)
if ~isnan(decay_log_slope)
    if     peak_latency_ms <= 40 && decay_log_slope <  -3.7;  classification = 'FS';
    elseif peak_latency_ms <  40 && decay_log_slope >  -3.7;  classification = 'BT';
    elseif peak_latency_ms >  40 && decay_log_slope >= -250;  classification = 'RF';
    else;                                                      classification = 'UNGROUPED';
    end
else
    classification = 'UNGROUPED';
end

ac_params.autocorr        = ac;
ac_params.lags_ms         = lags_s * 1000;
ac_params.peak_latency_ms = peak_latency_ms;
ac_params.decay_log_slope = decay_log_slope;
ac_params.classification  = classification;

if plot_it
    figure('Position', [50 50 1100 800]);
    subplot(2, 1, 1);
    plot(lags_s * 1000, ac, 'k', 'LineWidth', 1.5); hold on;
    plot(peak_latency_ms, ac(peak_idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    if ~isnan(decay_log_slope)
        plot(lags_s(fit_range) * 1000, ac(fit_range), 'b-', 'LineWidth', 2);
    end
    xline(40, '--k', 'LineWidth', 1.5, 'Alpha', 0.5, 'HandleVisibility', 'off');
    xlabel('Lag (ms)'); ylabel('Normalized Autocorrelation');
    title(sprintf('Peak Latency: %.1f ms, Decay log(slope): %.2f, Class: %s', ...
        peak_latency_ms, decay_log_slope, classification));
    xlim([0 max_lag_ms]); grid on; box off;

    subplot(2, 1, 2);
    valid_ac = ac > 0;
    plot(lags_s(valid_ac) * 1000, log(ac(valid_ac)), 'k', 'LineWidth', 1.5); hold on;
    if ~isnan(decay_log_slope)
        fit_range_valid = fit_range(ac(fit_range) > 0);
        if ~isempty(fit_range_valid)
            plot(lags_s(fit_range_valid) * 1000, log(ac(fit_range_valid)), 'b-', 'LineWidth', 2);
        end
    end
    xlabel('Lag (ms)'); ylabel('log(Autocorrelation)');
    title('Log Scale (for decay visualization)');
    xlim([0 max_lag_ms]); grid on; box off;
    pause(0.5); close;
end
end
