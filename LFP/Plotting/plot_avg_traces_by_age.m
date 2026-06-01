function plot_avg_traces_by_age(SR_DATA, T_VEC, T_PRM, PARAMS, varargin)
%  Mean +/- SEM LFP traces per Region x Window x
% Layer-variant. G
% SS 2026

PRM_amplitude_var = 'lowest';
PRM_x_mono_ms     = [-5, 60];
PRM_x_poly_ms     = [0, 100];
PRM_old_darken    = 0.45;
PRM_line_width    = 1.5;
Extract_varargin;

trace_field = sprintf('trace_per_stim_%s', PRM_amplitude_var);
out_dir     = fullfile(T_PRM.FIG_DIR, 'AvgTraces', PRM_amplitude_var);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

plot_specs = {
    'PL', 'mono', 'Mono', 'L5';
    'PL', 'poly', 'Poly', 'L2';
    'PL', 'poly', 'Poly', 'L5';
    'PL', 'poly', 'Poly', 'Combined';
    'IL', 'mono', 'Mono', 'L2';
    'IL', 'mono', 'Mono', 'L5';
    'IL', 'mono', 'Mono', 'Combined';
    'IL', 'poly', 'Poly', 'L2';
    'IL', 'poly', 'Poly', 'L5';
    'IL', 'poly', 'Poly', 'Combined';
};

for i = 1:size(plot_specs, 1)
    region  = plot_specs{i, 1};
    win_fld = plot_specs{i, 2};
    win_lbl = plot_specs{i, 3};
    variant = plot_specs{i, 4};
    if strcmp(win_fld, 'mono'); x_range = PRM_x_mono_ms;
    else;                       x_range = PRM_x_poly_ms;
    end
    plot_one_figure(SR_DATA, T_VEC, T_PRM, PARAMS, region, win_fld, win_lbl, ...
        variant, trace_field, x_range, PRM_old_darken, PRM_line_width, out_dir);
end
fprintf('  Avg-trace figures saved to: %s\n', out_dir);
end


function plot_one_figure(SR_DATA, T_VEC, T_PRM, PARAMS, region, win_fld, win_lbl, ...
    variant, trace_field, x_range, old_darken, line_width, out_dir)
base_col = PARAMS.colors.(region);
col_Y    = base_col;
col_O    = base_col * (1 - old_darken);

iHC = collect_per_rat_traces(SR_DATA, T_PRM, region, win_fld, variant, trace_field, 'iHC');
vHC = collect_per_rat_traces(SR_DATA, T_PRM, region, win_fld, variant, trace_field, 'vHC');
pool = nanmean(cat(3, iHC, vHC), 3);

ageY = strcmp(T_PRM.AGES, 'Y');
ageO = strcmp(T_PRM.AGES, 'O');

panels = {
    iHC(ageY, :),  iHC(ageO, :),  'iHC';
    vHC(ageY, :),  vHC(ageO, :),  'vHC';
    pool(ageY, :), pool(ageO, :), 'iHC + vHC';
};

fig = figure('Position', [100 100 1200 380], 'Color', 'w', 'Visible', 'off');
for s = 1:3
    subplot(1, 3, s); hold on;
    plot_mean_sem(T_VEC, panels{s, 1}, col_Y, line_width);
    plot_mean_sem(T_VEC, panels{s, 2}, col_O, line_width);
    yline(0, ':', 'Color', [0.6 0.6 0.6]);
    xlim(x_range); xlabel('Time (ms)'); ylabel('LFP (mV)');
    title(panels{s, 3}); box off;
end
sgtitle(sprintf('%s %s - %s   (Y light, O dark; mean \\pm SEM)', region, win_lbl, variant));
saveas(fig, fullfile(out_dir, sprintf('AvgTrace_%s_%s_%s.png', region, win_lbl, variant)));
close(fig);
end


function M = collect_per_rat_traces(SR_DATA, T_PRM, region, win_fld, variant, trace_field, sess_name)
% [n_rats x nT]: per-rat = mean across valid stims of per-stim channel-set
% trace. Rats with no data -> NaN row.
n_rats = numel(T_PRM.RATS);
nT = NaN; M = [];

for r = 1:n_rats
    ratFld = sprintf('Rat%s', T_PRM.RATS{r});
    if ~isfield(SR_DATA, ratFld); continue; end
    trace_r = get_variant_trace(SR_DATA.(ratFld), region, win_fld, variant, trace_field, sess_name);
    if isempty(trace_r); continue; end
    if isnan(nT); nT = numel(trace_r); M = nan(n_rats, nT); end
    if numel(trace_r) ~= nT; continue; end
    M(r, :) = trace_r(:)';
end
if isempty(M); M = nan(n_rats, fallback_nT(SR_DATA)); end
end


function trace = get_variant_trace(rat_entry, region, win_fld, variant, trace_field, sess_name)
% [1 x nT] per rat. Combined = per-rat mean of L2 and L5.
trace = [];
L2_sh = NaN; L5_sh = NaN;
if isfield(rat_entry, 'L2_shank'); L2_sh = rat_entry.L2_shank; end
if isfield(rat_entry, 'L5_shank'); L5_sh = rat_entry.L5_shank; end

switch variant
    case 'L2'; trace = pull_one_shank_trace(rat_entry, L2_sh, region, win_fld, trace_field, sess_name);
    case 'L5'; trace = pull_one_shank_trace(rat_entry, L5_sh, region, win_fld, trace_field, sess_name);
    case 'Combined'
        t2 = pull_one_shank_trace(rat_entry, L2_sh, region, win_fld, trace_field, sess_name);
        t5 = pull_one_shank_trace(rat_entry, L5_sh, region, win_fld, trace_field, sess_name);
        if     ~isempty(t2) && ~isempty(t5) && numel(t2) == numel(t5)
            trace = mean([t2(:)'; t5(:)'], 1, 'omitnan');
        elseif ~isempty(t2); trace = t2;
        elseif ~isempty(t5); trace = t5;
        end
end
end


function trace = pull_one_shank_trace(rat_entry, shIdx, region, win_fld, trace_field, sess_name)
% Per-rat = mean across valid stims of per-stim channel-set-averaged trace.
trace = [];
if isnan(shIdx); return; end
shFld = sprintf('Shank%d', shIdx);
if ~isfield(rat_entry, shFld) || ~isfield(rat_entry.(shFld), 'sessions'); return; end

sess_struct = rat_entry.(shFld).sessions;
sess_idx = NaN;
for s = 1:numel(sess_struct)
    if isfield(sess_struct(s), 'name') && strcmp(sess_struct(s).name, sess_name)
        sess_idx = s; break;
    end
end
if isnan(sess_idx); return; end

if ~isfield(sess_struct(sess_idx), region); return; end
reg = sess_struct(sess_idx).(region);
if ~isfield(reg, win_fld); return; end
pk = reg.(win_fld);
if isempty(pk) || ~isfield(pk, 'peak_detected') || ~pk.peak_detected; return; end
if ~isfield(pk, trace_field); return; end

ts = pk.(trace_field);
if isempty(ts); return; end
ok = ~all(isnan(ts), 2);
if ~any(ok); return; end
trace = mean(ts(ok, :), 1, 'omitnan');
end


function plot_mean_sem(t_ms, M, col, lw)
if isempty(M); return; end
ok = ~all(isnan(M), 2);
M  = M(ok, :);
if isempty(M); return; end
n   = size(M, 1);
mu  = mean(M, 1, 'omitnan');
sem = std(M, 0, 1, 'omitnan') ./ sqrt(n);
fill([t_ms(:)', fliplr(t_ms(:)')], [mu + sem, fliplr(mu - sem)], col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_ms, mu, 'Color', col, 'LineWidth', lw);
end


function nT = fallback_nT(SR_DATA)
% nT used when no rat had data; ensures NaN matrix has correct width.
rats = fieldnames(SR_DATA);
rats = rats(startsWith(rats, 'Rat'));
nT = 1;
for r = 1:numel(rats)
    rat_entry = SR_DATA.(rats{r});
    shanks = fieldnames(rat_entry);
    shanks = shanks(startsWith(shanks, 'Shank'));
    for s = 1:numel(shanks)
        if ~isfield(rat_entry.(shanks{s}), 'sessions'); continue; end
        ss = rat_entry.(shanks{s}).sessions;
        for k = 1:numel(ss)
            if isfield(ss(k), 'PL') && isfield(ss(k).PL, 'mono') ...
                    && isfield(ss(k).PL.mono, 'best_channel_traces') ...
                    && ~isempty(ss(k).PL.mono.best_channel_traces)
                nT = size(ss(k).PL.mono.best_channel_traces, 2); return;
            end
        end
    end
end
end
