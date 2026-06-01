function plot_mono_verification(MONO_DATA, T_VEC, PARAMS, FIG_DIR)


if ~exist(FIG_DIR, 'dir'); mkdir(FIG_DIR); end

[~, p_lo] = min(abs(T_VEC - (-5)));
[~, p_hi] = min(abs(T_VEC - 80));
t_plot    = T_VEC(p_lo:p_hi);

c_IL = PARAMS.colors.IL;

rats = fieldnames(MONO_DATA);
rats = rats(startsWith(rats, 'Rat'));

for iRat = 1:numel(rats)
    ratFld = rats{iRat};
    rat_id = MONO_DATA.(ratFld).rat_id;
    age    = MONO_DATA.(ratFld).age;
    L2_sh  = MONO_DATA.(ratFld).L2_shank;
    L5_sh  = MONO_DATA.(ratFld).L5_shank;

    for sess = 1:2
        sess_name = resolve_session_name(MONO_DATA.(ratFld), sess);
        if isempty(sess_name); continue; end

        fig = figure('Position', [50, 50, 1800, 450], 'Visible', 'off');
        tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

        for shIdx = 1:4
            shFld = sprintf('Shank%d', shIdx);
            nexttile(shIdx);
            if ~isfield(MONO_DATA.(ratFld), shFld); continue; end
            if ~isfield(MONO_DATA.(ratFld).(shFld), 'sessions'); continue; end
            sessions = MONO_DATA.(ratFld).(shFld).sessions;
            if numel(sessions) < sess; continue; end
            sess_data = sessions(sess);

            layer_tag = build_layer_tag(shIdx, L2_sh, L5_sh);
            plot_one_panel(sess_data.IL.mono, sess_data.valid_stims, ...
                sess_data.IL.strong_stims, sess_data.stim_currents, ...
                sess_data.IL.channel_mode, sess_data.mono_window, ...
                sess_data.window_override, sess_data.IL.frac_walled, ...
                t_plot, p_lo, p_hi, c_IL, sprintf('Sh%d%s', shIdx, layer_tag));
        end

        sgtitle(sprintf('Rat %s (%s) - %s', rat_id, age, sess_name), 'FontSize', 13);
        saveas(fig, fullfile(FIG_DIR, sprintf('VerifyMono_Rat%s_%s.png', rat_id, sess_name)));
        close(fig);
    end
    fprintf('  Plotted Rat %s\n', rat_id);
end
end


function plot_one_panel(pk, valid_stims, strong_stims, stim_currents, ...
    channel_mode, mono_win, has_override, frac_walled, ...
    t_plot, p_lo, p_hi, base_color, ttl)
hold on;
yline(0, 'k:', 'LineWidth', 0.5);
xline(mono_win(1), '--', 'Color', base_color * 0.6, 'LineWidth', 0.7);
xline(mono_win(2), '--', 'Color', base_color * 0.6, 'LineWidth', 0.7);

if isempty(pk) || ~isstruct(pk) || ~isfield(pk, 'peak_detected') || ~pk.peak_detected
    title([ttl ' - IL (no peak)'], 'FontSize', 10);
    xlabel('Time (ms)', 'FontSize', 9); ylabel('Voltage (mV)', 'FontSize', 9);
    xlim([t_plot(1), t_plot(end)]); set(gca, 'FontSize', 9);
    return;
end
if isempty(pk.best_channel_traces)
    title([ttl ' - IL (no traces)'], 'FontSize', 10);
    xlim([t_plot(1), t_plot(end)]); return;
end

traces    = pk.best_channel_traces(:, p_lo:p_hi);
n_stim    = size(traces, 1);
cmap      = make_stim_gradient(base_color, n_stim);
is_strong = ismember(valid_stims(:), strong_stims(:));

for i = 1:n_stim
    if is_strong(i); continue; end
    plot(t_plot, traces(i, :), '-', 'Color', cmap(i, :), 'LineWidth', 0.6);
end
for i = 1:n_stim
    if ~is_strong(i); continue; end
    plot(t_plot, traces(i, :), '-', 'Color', cmap(i, :), 'LineWidth', 1.6);
end

plot(pk.best_channel_latency, pk.best_channel_amplitude, 'o', ...
    'MarkerFaceColor', base_color, 'MarkerEdgeColor', base_color * 0.4, 'MarkerSize', 4);

strong_uA = round(stim_currents(strong_stims(:)));
ov_tag    = ''; if has_override; ov_tag = '*'; end
n_excl    = numel(pk.excluded_channels);
excl_tag  = ''; if n_excl > 0; excl_tag = sprintf(' ex=%d', n_excl); end
wall_tag  = ''; if ~isnan(frac_walled); wall_tag = sprintf(' wall=%.2f', frac_walled); end

title(sprintf('%s%s ch%d [%s] win=[%g %g]%s%s\nstrong: [%s]uA', ...
    ttl, ov_tag, pk.best_channel, channel_mode, mono_win(1), mono_win(2), ...
    wall_tag, excl_tag, num2str(strong_uA(:)')), 'FontSize', 9);

xlabel('Time (ms)', 'FontSize', 9); ylabel('Voltage (mV)', 'FontSize', 9);
xlim([t_plot(1), t_plot(end)]); set(gca, 'FontSize', 9);
end


function tag = build_layer_tag(shIdx, L2_sh, L5_sh)
if     ~isnan(L2_sh) && shIdx == L2_sh; tag = '[L2]';
elseif ~isnan(L5_sh) && shIdx == L5_sh; tag = '[L5]';
else;                                    tag = '';
end
end


function nm = resolve_session_name(ratData, sess)
nm = '';
shanks = fieldnames(ratData);
shanks = shanks(startsWith(shanks, 'Shank'));
for s = 1:numel(shanks)
    if ~isfield(ratData.(shanks{s}), 'sessions'); continue; end
    sd = ratData.(shanks{s}).sessions;
    if numel(sd) >= sess && ~isempty(sd(sess).name); nm = sd(sess).name; return; end
end
end
