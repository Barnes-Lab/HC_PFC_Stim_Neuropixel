function plot_slope_verification(SR_DATA, T_VEC, PARAMS, FIG_DIR)
% Verified detection of peak of response slope Pe
% PL-mono, PL-poly, IL-mono, IL-poly). 
% best-channel traces light todark with increasing stim intensity i
% SS 2026

if ~exist(FIG_DIR, 'dir'); mkdir(FIG_DIR); end

[~, p_lo] = min(abs(T_VEC - (-5)));
[~, p_hi] = min(abs(T_VEC - 150));
t_plot    = T_VEC(p_lo:p_hi);

c_PL = PARAMS.colors.PL;
c_IL = PARAMS.colors.IL;

rats = fieldnames(SR_DATA);
rats = rats(startsWith(rats, 'Rat'));

for iRat = 1:numel(rats)
    ratFld = rats{iRat};
    rat_id = SR_DATA.(ratFld).rat_id;
    age    = SR_DATA.(ratFld).age;
    L2_sh  = SR_DATA.(ratFld).L2_shank;
    L5_sh  = SR_DATA.(ratFld).L5_shank;

    layer_shanks = [L2_sh, L5_sh];
    layer_tags   = {'L2', 'L5'};

    sess_names = resolve_session_names(SR_DATA.(ratFld), L2_sh);
    for sess = 1:numel(sess_names)
        sess_name = sess_names{sess};
        if isempty(sess_name); continue; end

        fig = figure('Position', [50, 50, 1800, 700], 'Visible', 'off');
        tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

        for li = 1:2
            shIdx = layer_shanks(li);
            shFld = sprintf('Shank%d', shIdx);
            if ~isfield(SR_DATA.(ratFld), shFld); continue; end
            sd = SR_DATA.(ratFld).(shFld).sessions(sess);
            base_idx = (li - 1) * 4;

            nexttile(base_idx + 1);
            plot_one_panel(sd.PL.mono, sd.PL.mono_window, t_plot, p_lo, p_hi, ...
                c_PL,       sprintf('Sh%d %s - PL mono', shIdx, layer_tags{li}));
            nexttile(base_idx + 2);
            plot_one_panel(sd.PL.poly, sd.PL.poly_window, t_plot, p_lo, p_hi, ...
                c_PL * 0.6, sprintf('Sh%d %s - PL poly', shIdx, layer_tags{li}));
            nexttile(base_idx + 3);
            plot_one_panel(sd.IL.mono, sd.IL.mono_window, t_plot, p_lo, p_hi, ...
                c_IL,       sprintf('Sh%d %s - IL mono', shIdx, layer_tags{li}));
            nexttile(base_idx + 4);
            plot_one_panel(sd.IL.poly, sd.IL.poly_window, t_plot, p_lo, p_hi, ...
                c_IL * 0.6, sprintf('Sh%d %s - IL poly', shIdx, layer_tags{li}));
        end

        sgtitle(sprintf('Rat %s (%s) - %s', rat_id, age, sess_name), 'FontSize', 13);
        saveas(fig, fullfile(FIG_DIR, sprintf('VerifySlope_Rat%s_%s.png', rat_id, sess_name)));
        close(fig);
    end
end
end


function plot_one_panel(pk, win, t_plot, p_lo, p_hi, base_color, ttl)
hold on;
yline(0, 'k:', 'LineWidth', 0.5);
xline(win(1), '--', 'Color', base_color * 0.5, 'LineWidth', 0.7);
xline(win(2), '--', 'Color', base_color * 0.5, 'LineWidth', 0.7);

if isempty(pk) || ~isstruct(pk) || ~isfield(pk, 'peak_detected') ...
        || ~pk.peak_detected || isempty(pk.best_channel_traces)
    title([ttl ' (no peak)'], 'FontSize', 9);
    xlabel('Time (ms)', 'FontSize', 8); ylabel('mV', 'FontSize', 8);
    xlim([t_plot(1), t_plot(end)]); set(gca, 'FontSize', 8);
    return;
end

traces = pk.best_channel_traces(:, p_lo:p_hi);
n_stim = size(traces, 1);
cmap   = make_stim_gradient(base_color, n_stim);
for i = 1:n_stim
    plot(t_plot, traces(i, :), '-', 'Color', cmap(i, :), 'LineWidth', 0.9);
end

plot(pk.best_channel_latency, pk.best_channel_amplitude, 'o', ...
    'MarkerFaceColor', base_color, 'MarkerEdgeColor', base_color * 0.3, 'MarkerSize', 5);

title(sprintf('%s\nch%d w[%g %g] n50=%d', ttl, pk.best_channel, win(1), win(2), ...
    pk.n_region_channels_50pct), 'FontSize', 8);
xlabel('Time (ms)', 'FontSize', 8); ylabel('mV', 'FontSize', 8);
xlim([t_plot(1), t_plot(end)]); set(gca, 'FontSize', 8);
end


function names = resolve_session_names(ratData, ref_shank)
names = {};
shFld = sprintf('Shank%d', ref_shank);
if ~isfield(ratData, shFld) || ~isfield(ratData.(shFld), 'sessions'); return; end
sd = ratData.(shFld).sessions;
for s = 1:numel(sd)
    if isfield(sd(s), 'name'); names{end+1, 1} = sd(s).name; end
end
end
