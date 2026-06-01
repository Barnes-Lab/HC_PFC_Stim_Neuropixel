function plot_within_rat(data, sess_types, comparison_type, metric, ageCols, ...
    figDir, y_label, recruit_thresh)
% PLOT_WITHIN_RAT  Within-rat paired plot (PYR only, neuron-weighted).
% Layouts:
%   'mono_poly'  1x3  All | PL | IL                  900x500, pubify 16
%   'PL_IL'      3x2  (iHC, vHC, All) x (mono, poly) 600x900, pubify 14
%   'iHC_vHC'    2x2  (PL, IL) x (mono, poly)        600x600, pubify 14
% Saves .png and .fig.
% SS 2026

[titles, fAs, fBs, layout, xtick_labels, fig_size, pubify_size] = get_layout(comparison_type);

rats_Y = unique_rats_for_age(data, sess_types, 'PYR', 'Y');
rats_O = unique_rats_for_age(data, sess_types, 'PYR', 'O');

figure('Position', [100, 100, fig_size(1), fig_size(2)]);
for k = 1:numel(titles)
    fA = fAs{k}; fB = fBs{k};
    vY_a = get_per_rat_value(data, sess_types, 'Y', fA{1}, fA{2}, fA{3}, metric, recruit_thresh, rats_Y);
    vY_b = get_per_rat_value(data, sess_types, 'Y', fB{1}, fB{2}, fB{3}, metric, recruit_thresh, rats_Y);
    vO_a = get_per_rat_value(data, sess_types, 'O', fA{1}, fA{2}, fA{3}, metric, recruit_thresh, rats_O);
    vO_b = get_per_rat_value(data, sess_types, 'O', fB{1}, fB{2}, fB{3}, metric, recruit_thresh, rats_O);

    ax = subplot(layout(1), layout(2), k);
    plot_within_rat_panel(ax, vY_a, vY_b, vO_a, vO_b, ageCols, ...
        xtick_labels, y_label, titles{k}, 'PRM_PUBIFY_SIZE', pubify_size);
end
sgtitle(sprintf('Within-Rat %s: %s [PYR]', comparison_label(comparison_type), metric), ...
    'Interpreter', 'none');

base = fullfile(figDir, sprintf('within_rat_%s_%s_PYR', comparison_type, metric));
saveas(gcf, [base '.png']);
saveas(gcf, [base '.fig']);
close(gcf);
end


function lbl = comparison_label(comparison_type)
switch comparison_type
    case 'mono_poly'; lbl = 'Mono vs Poly';
    case 'PL_IL';     lbl = 'PL vs IL';
    case 'iHC_vHC';   lbl = 'iHC vs vHC';
    otherwise;        lbl = comparison_type;
end
end


function [titles, fAs, fBs, layout, xtick_labels, fig_size, pubify_size] = get_layout(comparison_type)
switch comparison_type
    case 'mono_poly'
        layout       = [1, 3];
        fig_size     = [900, 500];
        pubify_size  = 16;
        xtick_labels = {'Mono', 'Poly'};
        titles = {'All', 'PL', 'IL'};
        fAs    = {{'all','avg','mono'}, {'all','PL','mono'}, {'all','IL','mono'}};
        fBs    = {{'all','avg','poly'}, {'all','PL','poly'}, {'all','IL','poly'}};

    case 'PL_IL'
        layout       = [3, 2];
        fig_size     = [600, 900];
        pubify_size  = 14;
        xtick_labels = {'PL', 'IL'};
        sess_levels  = {'iHC','vHC','all'};
        sess_disp    = {'iHC','vHC','All'};
        wins         = {'mono','poly'};
        n = 6;
        titles = cell(1, n); fAs = cell(1, n); fBs = cell(1, n);
        k = 0;
        for si = 1:3
            for wi = 1:2
                k = k + 1;
                titles{k} = sprintf('%s - %s', sess_disp{si}, wins{wi});
                fAs{k}    = {sess_levels{si}, 'PL', wins{wi}};
                fBs{k}    = {sess_levels{si}, 'IL', wins{wi}};
            end
        end

    case 'iHC_vHC'
        layout       = [2, 2];
        fig_size     = [600, 600];
        pubify_size  = 14;
        xtick_labels = {'iHC', 'vHC'};
        regions      = {'PL', 'IL'};
        wins         = {'mono', 'poly'};
        n = 4;
        titles = cell(1, n); fAs = cell(1, n); fBs = cell(1, n);
        k = 0;
        for ri = 1:2
            for wi = 1:2
                k = k + 1;
                titles{k} = sprintf('%s - %s', regions{ri}, wins{wi});
                fAs{k}    = {'iHC', regions{ri}, wins{wi}};
                fBs{k}    = {'vHC', regions{ri}, wins{wi}};
            end
        end

    otherwise
        error('Unknown comparison_type: %s', comparison_type);
end
end
