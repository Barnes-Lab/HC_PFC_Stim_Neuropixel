function plot_i80_verification(rat_id, sess_name, shank_results, ...
                                chosen_shank, I80_final, save_path, varargin)
% 1x4 verification [plots showing the IO curve and fit on each shank
% for each rat per sess
% 180 stim is marked along with chosen shank and R2 value of fit
% SS 2026

%% PARAMETERS
PRM_chosen_source = '';
Extract_varargin;

orange = [0.88, 0.35, 0.11];
blue   = [0.00, 0.40, 0.80];
grey   = [0.55, 0.55, 0.55];

flat_sources   = {'manual_override_flat','forced_flat','forced_flat_other_shank'};
chosen_is_flat = any(strcmp(PRM_chosen_source, flat_sources));

logistic = @(Amin, Amax, k, I50, xv) ...
    Amin + (Amax - Amin) ./ (1 + exp(-k .* (xv - I50)));
expsat   = @(Amax, k, x0, xv) ...
    Amax .* (1 - exp(-k .* (xv - x0)));

fig = figure('Position', [100, 100, 1800, 420], 'Visible', 'off');

for s = 1:4
    subplot(1, 4, s); hold on;
    sr = shank_results(s);
    layer_str = sr.layer; if isempty(layer_str); layer_str = '-'; end

    if ~sr.has_data
        text(0.5, 0.5, sprintf('no data\n(%s)', sr.Status), ...
            'Units','normalized','HorizontalAlignment','center', ...
            'Interpreter','none','Color',grey,'FontSize',11);
        title(sprintf('Sh%d %s', s, layer_str));
        xlim([50 650]); pubify_figure_axis_robust(14, 14);
        format_chosen(s, chosen_shank, blue);
        continue;
    end

    plot(sr.x, sr.y, 'o', 'Color', grey, 'MarkerFaceColor', grey, ...
        'MarkerSize', 6);

    % Fit curve
    xs = linspace(min(sr.x), max(sr.x), 100);
    f  = sr.fit;
    if (s == chosen_shank) && chosen_is_flat
        yflat = mean(sr.y);
        plot([min(sr.x), max(sr.x)], [yflat, yflat], '-', ...
            'Color', orange, 'LineWidth', 2);
    elseif contains(f.Method, 'logistic')
        plot(xs, logistic(f.Amin, f.Amax, f.k, f.I50, xs), '-', ...
            'Color', orange, 'LineWidth', 2);
    elseif contains(f.Method, 'exp_sat')
        plot(xs, expsat(f.Amax, f.k, f.x0, xs), '-', ...
            'Color', orange, 'LineWidth', 2);
    end

    % I80 line
    if s == chosen_shank
        I80_show = I80_final;
    else
        I80_show = f.I80;
    end
    if ~isnan(I80_show) && I80_show > 0
        yL = ylim;
        plot([I80_show, I80_show], [yL(1), yL(2)], '--', ...
            'Color', 'r', 'LineWidth', 1.5);
        text(I80_show, yL(1) + 0.06 * (yL(2) - yL(1)), ...
            sprintf('%.0f', I80_show), 'Color', 'r', 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end

    title(sprintf('Sh%d %s | %s | R^2=%.2f', s, layer_str, ...
        strrep(f.Method, '_', ' '), f.R2));
    xlabel('Current (\muA)'); ylabel('Amp (mV)');
    xlim([50 650]); box off;
    pubify_figure_axis_robust(14, 14);
    format_chosen(s, chosen_shank, blue);
end

if isnan(chosen_shank)
    src_str = 'NO SHANK CHOSEN';
else
    src_str = sprintf('chosen Sh%d (I80=%.0f \\muA)', chosen_shank, I80_final);
end
sgtitle(sprintf('Rat %s — %s — %s', rat_id, sess_name, src_str), ...
    'FontSize', 14, 'FontWeight', 'bold');

saveas(fig, save_path);
close(fig);

end

% Inline format easier
function format_chosen(s, chosen_shank, blue)
if s == chosen_shank
    ax = gca;
    ax.LineWidth      = 3;
    ax.XColor         = blue;
    ax.YColor         = blue;
    ax.Title.Color    = blue;
    ax.Title.FontWeight = 'bold';
end
end
