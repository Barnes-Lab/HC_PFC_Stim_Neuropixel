function create_all_I80_figures(data, COLORS, T_PRM)
% CREATE_ALL_I80_FIGURES - Generate all publication figures using I80-selected data
%
% Uses existing plotting functions from create_all_fEPSP_figures
% For Combined analysis, adds Figure 4: Within-rat L2 vs L5 comparison

fprintf('=== Creating Figures (I80-Based Stimulus Selection) ===\n');

% Detect if this is Combined analysis
is_combined = isfield(data, 'layers') && ...
              any(strcmp(data.layers, 'L2')) && ...
              any(strcmp(data.layers, 'L5'));

%% 1: Within-Rat Mono vs Poly
fprintf('Creating Figure 1: Mono vs Poly...\n');
fig1 = figure('Position', [100 100 600 500]);
plot_mono_vs_poly_with_means(data, COLORS);
saveas(fig1, fullfile(T_PRM.PUB_FIG_DIR, 'Figure1_Mono_vs_Poly.png'));
saveas(fig1, fullfile(T_PRM.PUB_FIG_DIR, 'Figure1_Mono_vs_Poly.fig'));
close(fig1);

%%  2: Within-Rat PL vs IL (2figs-amplitude and latency)
fprintf('Creating Figure 2: PL vs IL...\n');
figs_before = findobj('Type', 'figure');
plot_PL_vs_IL_with_means(data, COLORS);
figs_after = findobj('Type', 'figure');
new_figs = setdiff(figs_after, figs_before);

for i = 1:length(new_figs)
    fig_title = get(new_figs(i), 'Name');
    if isempty(fig_title)
        children = get(new_figs(i), 'Children');
        for j = 1:length(children)
            if isprop(children(j), 'Title')
                title_obj = get(children(j), 'Title');
                if ~isempty(title_obj)
                    fig_title = get(title_obj, 'String');
                    break;
                end
            end
        end
    end
    
    if contains(fig_title, 'Amplitude')
        saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure2A_PL_vs_IL_Amplitude.png'));
        saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure2A_PL_vs_IL_Amplitude.fig'));
    else
        saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure2B_PL_vs_IL_Latency.png'));
        saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure2B_PL_vs_IL_Latency.fig'));
    end
    close(new_figs(i));
end

%% 3: Within-Rat iHC vs vHC 
fprintf('Creating Figure 3: iHC vs vHC...\n');
figs_before = findobj('Type', 'figure');
plot_iHC_vs_vHC_with_means(data, COLORS);
figs_after = findobj('Type', 'figure');
new_figs = setdiff(figs_after, figs_before);
for i = 1:length(new_figs)
    fig_title = get(new_figs(i), 'Name');
    if isempty(fig_title)
        children = get(new_figs(i), 'Children');
        for j = 1:length(children)
            if isprop(children(j), 'Title')
                title_obj = get(children(j), 'Title');
                if ~isempty(title_obj)
                    fig_title = get(title_obj, 'String');
                    break;
                end
            end
        end
    end
    if contains(fig_title, 'Amplitude')
        saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure3A_iHC_vs_vHC_Amplitude.png'));
        saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure3A_iHC_vs_vHC_Amplitude.fig'));
    else
        saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure3B_iHC_vs_vHC_Latency.png'));
        saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure3B_iHC_vs_vHC_Latency.fig'));
    end
    close(new_figs(i));
end

%% 4: Within-Rat L2 vs L5 
if is_combined
    fprintf('Creating Figure 4: L2 vs L5 (Combined analysis)...\n');
    figs_before = findobj('Type', 'figure');
    plot_L2_vs_L5_with_means(data, COLORS);
    figs_after = findobj('Type', 'figure');
    new_figs = setdiff(figs_after, figs_before);
    for i = 1:length(new_figs)
        fig_title = get(new_figs(i), 'Name');
        if isempty(fig_title)
            children = get(new_figs(i), 'Children');
            for j = 1:length(children)
                if isprop(children(j), 'Title')
                    title_obj = get(children(j), 'Title');
                    if ~isempty(title_obj)
                        fig_title = get(title_obj, 'String');
                        break;
                    end
                end
            end
        end
        if contains(fig_title, 'Amplitude')
            saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure4A_L2_vs_L5_Amplitude.png'));
            saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure4A_L2_vs_L5_Amplitude.fig'));
        else
            saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure4B_L2_vs_L5_Latency.png'));
            saveas(new_figs(i), fullfile(T_PRM.PUB_FIG_DIR, 'Figure4B_L2_vs_L5_Latency.fig'));
        end
        close(new_figs(i));
    end
end

%% Bar plots old
fprintf('Creating grouped bar plots...\n');
plot_bars_grouped_age(data, COLORS, T_PRM.PUB_FIG_DIR);

fprintf('\n=== All figures complete! ===\n');

end
