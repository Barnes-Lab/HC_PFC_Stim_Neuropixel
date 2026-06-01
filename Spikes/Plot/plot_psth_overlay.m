function plot_psth_overlay(data, bin_centers, xlim_ms, neuron_types, ageCols, ...
    mode, y_label, sg_title, fname_no_ext)
% PLOT_PSTH_OVERLAY  2-row figure (one row per neuron type), 4 traces per
% panel (Y-PL, Y-IL, O-PL, O-IL) overlaid.
%
% Modes:
%   'pooled'      uses res.pooled_bsl, n = n_neurons (biased by N)
%   'perRat_bsl'  uses res.rat_bsl,    n = n_rats   (per-rat dFR, unbiased)
% SS 2026

ages    = {'Y', 'O'};
regions = {'PL', 'IL'};
ls      = struct('PL', '-', 'IL', '--');

figure('Position', [100, 100, 1000, 1000]);
for ti = 1:length(neuron_types)
    nt = neuron_types{ti};
    subplot(2, 1, ti); hold on;

    for ai = 1:2
        a   = ages{ai};
        col = ageCols.(a);
        for ri = 1:2
            r   = regions{ri};
            res = data.(nt).(a).(r);
            switch mode
                case 'pooled'
                    M = res.pooled_bsl; n_lbl = sprintf('n=%d', res.n_neurons);
                case 'perRat_bsl'
                    M = res.rat_bsl;    n_lbl = sprintf('n_{rats}=%d', res.n_rats);
                otherwise
                    error('Unknown mode: %s (use ''pooled'' or ''perRat_bsl'')', mode);
            end
            if size(M, 1) == 0; continue; end

            mp = mean(M, 1, 'omitnan');
            sp = std(M, 0, 1, 'omitnan') / sqrt(size(M, 1));
            plot(bin_centers, mp, 'LineStyle', ls.(r), 'Color', col, 'LineWidth', 3, ...
                'DisplayName', sprintf('%s-%s (%s)', a, r, n_lbl));
            fill([bin_centers, fliplr(bin_centers)], ...
                [mp + sp, fliplr(mp - sp)], col, 'FaceAlpha', 0.2, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end

    xline(0,   'Color', '#A2142F', 'LineWidth', 3, 'HandleVisibility', 'off');
    xline(35,  'k:', 'LineWidth', 2, 'HandleVisibility', 'off');
    xline(100, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');
    xlim(xlim_ms);
    xlabel('Time from stim (ms)');
    ylabel(y_label);
    title(sprintf('%s Neurons', nt));
    legend('Location', 'bestoutside');
    box off;
    if exist('pubify_figure_axis_robust', 'file'); pubify_figure_axis_robust(20, 20); end
end

sgtitle(sg_title);
saveas(gcf, [fname_no_ext '.png']);
saveas(gcf, [fname_no_ext '.fig']);
close(gcf);
end
