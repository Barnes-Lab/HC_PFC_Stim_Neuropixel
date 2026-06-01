function plot_i80_paper_perrat(I80_DATA, T_PRM, output_dir)
% SS 2026


if ~exist(output_dir, 'dir'); mkdir(output_dir); end

orange = [0.95, 0.45, 0.10];
red    = [0.85, 0.10, 0.10];

logistic = @(Amin, Amax, k, I50, xv) ...
    Amin + (Amax - Amin) ./ (1 + exp(-k .* (xv - I50)));
expsat   = @(Amax, k, x0, xv) Amax .* (1 - exp(-k .* (xv - x0)));

for iRat = 1:numel(T_PRM.RATS)
    rat_id = T_PRM.RATS{iRat};
    age    = T_PRM.AGES{iRat};
    ratFld = sprintf('Rat%s', rat_id);
    if ~isfield(I80_DATA, ratFld); continue; end

    fig = figure('Position', [100, 100, 1100, 480], 'Visible', 'off');

    sites = {'iHC', 'vHC'};
    for s = 1:2
        subplot(1, 2, s); hold on;
        site = sites{s};

        if ~isfield(I80_DATA.(ratFld), site)
            text(0.5, 0.5, 'no session', 'Units','normalized', ...
                'HorizontalAlignment','center','FontSize',14);
            title(site, 'FontWeight', 'bold');
            xlim([50 650]); pubify_figure_axis_robust(18, 18);
            continue;
        end

        d  = I80_DATA.(ratFld).(site);
        sh = d.chosen_shank;
        if isnan(sh) || sh < 1 || sh > numel(d.shank_results); continue; end
        sr = d.shank_results(sh);

        % --- scatter
        if sr.has_data
            plot(sr.x, sr.y, 'o', 'Color', orange, ...
                'MarkerFaceColor', orange, 'MarkerSize', 9, ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end

        % --- best fit curve (clamped to data x range to avoid extrapolation blow-up)
        if sr.has_data
            xs = linspace(min(sr.x), max(sr.x), 200);
            f  = sr.fit;
            if contains(f.Method, 'logistic')
                plot(xs, logistic(f.Amin, f.Amax, f.k, f.I50, xs), '-', ...
                    'Color', orange, 'LineWidth', 3);
            elseif contains(f.Method, 'exp_sat')
                plot(xs, expsat(f.Amax, f.k, f.x0, xs), '-', ...
                    'Color', orange, 'LineWidth', 3);
            elseif strcmp(d.chosen_source, 'manual_override_flat') || ...
                   strcmp(d.chosen_source, 'forced_flat') || ...
                   strcmp(d.chosen_source, 'forced_flat_other_shank')
                yflat = mean(sr.y);
                plot([min(sr.x) max(sr.x)], [yflat yflat], '-', ...
                    'Color', orange, 'LineWidth', 3);
            end
        end

        % Line
        I80v = d.I80;
        if ~isnan(I80v) && I80v > 0
            yL = ylim;
            plot([I80v, I80v], [yL(1), yL(2)], '--', ...
                'Color', red, 'LineWidth', 2);
            text(I80v, yL(1) + 0.05 * (yL(2) - yL(1)), ...
                sprintf(' I80 = %.0f \\muA', I80v), ...
                'Color', red, 'FontWeight', 'bold', ...
                'FontSize', 14, 'HorizontalAlignment', 'left');
        end

        title(site, 'FontWeight', 'bold');
        xlabel('Current (\muA)');
        ylabel('Amplitude (mV)');
        xlim([50 650]);
        pubify_figure_axis_robust(18, 18);
        box off;
    end

    sgtitle(sprintf('Rat %s (%s)', rat_id, age), ...
        'FontSize', 18, 'FontWeight', 'bold');

    saveas(fig, fullfile(output_dir, sprintf('Rat%s_I80_paper.png', rat_id)));
    close(fig);
end


end
