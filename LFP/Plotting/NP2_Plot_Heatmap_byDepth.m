function NP2_Plot_Heatmap_byDepth(T_PRM)
% Plots heatmap of average response across all stims
% Creates per-rat heatmaps and combined session averages for L2 and L5
fprintf('Plotting Average Heatmap across conditions')
gamma = 1;
plot_best_sh = true;


PRM_IMG_DIR = fullfile(T_PRM.IMG_DIR,"HeatMap_Average");
if ~exist(PRM_IMG_DIR, "dir"), mkdir(PRM_IMG_DIR); end

tbl = readtable(fullfile(T_PRM.STORAGE_DIR,'ALL_STIM_TIMES_KEEP_CHECK.csv'));

firstRat = T_PRM.RATS{1};
load(fullfile(T_PRM.STORAGE_DIR, firstRat, sprintf('%s_DATA.mat',firstRat)), ...
    'PRM','LF_Shank','STIM');
ts = -PRM.time_extract_preStim*1e3;
te = PRM.time_extract_postStim*1e3;
t_all = linspace(ts, te, size(LF_Shank(1).V{1},3));
[~, xs] = min(abs(t_all + 50));
[~, xe] = min(abs(t_all - 150));
t = t_all(xs:xe);
nT = numel(t);

%% Get shank map
shankMap = readtable(fullfile(T_PRM.STORAGE_DIR,'ShankInfo','ShankMap.csv'),...
    'VariableNamingRule','preserve');
shankTypes = {'L2','L5'};

%% Plot per-rat heatmaps for L2 and L5
maxChannels = 94;
nRats = numel(T_PRM.RATS);

for sType = 1:2
    shankLabel = shankTypes{sType};
    
    % Storage for combined sessions by age
    VsumY_combined = zeros(maxChannels, nT);
    VsumO_combined = zeros(maxChannels, nT);
    nY_combined = 0;
    nO_combined = 0;
    
    for sess = 1:2
        for r = 1:nRats
            rat = T_PRM.RATS{r};
            age = T_PRM.AGES{r};
            
            % Valid stims
            rows = tbl(strcmp(string(tbl.Rat), rat) & ...
                strcmp(tbl.Session, STIM.epochs{sess}), :);
            ex = false(1, numel(STIM.sorted_current_uA));
            for k = 1:numel(STIM.sorted_current_uA)
                idx = find(rows.Current_uA == round(STIM.sorted_current_uA(k)), 1);
                if ~isempty(idx) && strcmp(rows.Keep{idx}, 'N')
                    ex(k) = true;
                end
            end
            valid = find(~ex);
            
            % Get shank index
            rowMap = strcmp(string(shankMap.Rat), rat);
            if ~any(rowMap)
                error('Rat %s missing from ShankMap.csv', rat);
            end
            shank_idx = shankMap.(shankLabel)(rowMap);
            
            % Load data
            d = load(fullfile(T_PRM.STORAGE_DIR, rat, sprintf('%s_DATA.mat',rat)), 'LF_Shank');
            Vb = squeeze(d.LF_Shank(shank_idx).V{sess}(valid, 1:maxChannels, xs:xe));
            meanCh94 = squeeze(mean(Vb,1));
            
            depths = d.LF_Shank(shank_idx).SortedDepths(1:maxChannels);
            contrast_Vmat = real(meanCh94 .^ gamma);
            
            % Per-rat heatmap
            figure;
            imagesc(t, depths, contrast_Vmat);
            axis xy;
            set(gca,'YDir','reverse');
            colormap('jet'); colorbar;
            clim([-0.8 0.8]);
            hold on;
            line([0 0], [min(depths) max(depths)], 'Color','w','LineWidth',3);
            xlabel('Time (ms)');
            ylabel('Depth (\mum)');
            title(sprintf('%s - %s - %s', rat, STIM.epochs{sess}, shankLabel));
            pubify_figure_axis_robust(14,14);
            img_name = sprintf('%s_%s_%s_Heatmap', rat, STIM.epochs{sess}, shankLabel);
            saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
            saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
            close;
            
            % Accumulate for combined average
            if strcmp(age,'Y')
                VsumY_combined = VsumY_combined + meanCh94;
                if sess == 2
                    nY_combined = nY_combined + 1;
                end
            else
                VsumO_combined = VsumO_combined + meanCh94;
                if sess == 2
                    nO_combined = nO_combined + 1;
                end
            end
        end
    end
    
    % Create combined session averages (iHC + vHC) for Young and Old
    depths_avg = zeros(maxChannels,1);
    for sh = 1:T_PRM.shank_no
        depths_avg = depths_avg + LF_Shank(sh).SortedDepths(1:maxChannels);
    end
    depths_avg = depths_avg / T_PRM.shank_no;
    
    % Young combined
    VmatY_combined = VsumY_combined / (nY_combined * 2);  % Divide by num rats * 2 sessions
    CmatY_combined = real(VmatY_combined .^ gamma);
    
    figure;
    imagesc(t, depths_avg, CmatY_combined);
    axis xy;
    set(gca,'YDir','reverse');
    colormap('jet'); colorbar;
    clim([-0.75 0.75]);
    hold on;
    line([0 0], [min(depths_avg) max(depths_avg)], 'Color','w','LineWidth',3);
    xlabel('Time (ms)');
    ylabel('Depth (\mum)');
    title(sprintf('Young - Combined (iHC+vHC) - %s', shankLabel));
    pubify_figure_axis_robust(14,14);
    img_name = sprintf('Young_Combined_iHC_vHC_%s', shankLabel);
    saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
    saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
    close;
    
    % Old combined
    VmatO_combined = VsumO_combined / (nO_combined * 2);  % Divide by num rats * 2 sessions
    CmatO_combined = real(VmatO_combined .^ gamma);
    
    figure;
    imagesc(t, depths_avg, CmatO_combined);
    axis xy;
    set(gca,'YDir','reverse');
    colormap('jet'); colorbar;
    clim([-0.75 0.75]);
    hold on;
    line([0 0], [min(depths_avg) max(depths_avg)], 'Color','w','LineWidth',3);
    xlabel('Time (ms)');
    ylabel('Depth (\mum)');
    title(sprintf('Old - Combined (iHC+vHC) - %s', shankLabel));
    pubify_figure_axis_robust(14,14);
    img_name = sprintf('Old_Combined_iHC_vHC_%s', shankLabel);
    saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
    saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
    close;
end

%% Plot best shank 
if plot_best_sh
    maxChannels = 94;
    depths = zeros(maxChannels,1);
    for sh = 1:T_PRM.shank_no
        depths = depths + LF_Shank(sh).SortedDepths(1:maxChannels);
    end
    depths = depths / T_PRM.shank_no;

    youngIdx = find(strcmp(T_PRM.AGES,'Y'));
    oldIdx = find(strcmp(T_PRM.AGES,'O'));
    nRats = numel(T_PRM.RATS);

    for sess = 1:2
        Vsum = zeros(maxChannels, nT);
        VsumY = zeros(maxChannels, nT);
        VsumO = zeros(maxChannels, nT);

        % All rats
        for r = 1:nRats
            rat = T_PRM.RATS{r};
            rows = tbl(strcmp(string(tbl.Rat), rat) & ...
                strcmp(tbl.Session, STIM.epochs{sess}), :);
            ex = false(1, numel(STIM.sorted_current_uA));
            for k = 1:numel(STIM.sorted_current_uA)
                idx = find(rows.Current_uA == round(STIM.sorted_current_uA(k)), 1);
                if ~isempty(idx) && strcmp(rows.Keep{idx}, 'N')
                    ex(k) = true;
                end
            end
            valid = find(~ex);

            d = load(fullfile(T_PRM.STORAGE_DIR, rat, sprintf('%s_DATA.mat',rat)), 'LF_Shank');
            mags = zeros(1, T_PRM.shank_no);
            for sh = 1:T_PRM.shank_no
                Vsh = squeeze(d.LF_Shank(sh).V{sess}(valid, :, xs:xe));
                flat = reshape(Vsh,[],size(Vsh,3));
                mags(sh) = abs(min(mean(flat,1)));
            end
            [~, best_shank] = max(mags);

            Vb = squeeze(d.LF_Shank(best_shank).V{sess}(valid, 1:maxChannels, xs:xe));
            meanCh94 = squeeze(mean(Vb,1));
            Vsum = Vsum + meanCh94;
        end
        
        % Young rats
        for r = youngIdx
            rat = T_PRM.RATS{r};
            rows = tbl(strcmp(string(tbl.Rat), rat) & strcmp(tbl.Session,STIM.epochs{sess}), :);
            ex = false(1, numel(STIM.sorted_current_uA));
            for k = 1:numel(STIM.sorted_current_uA)
                idx = find(rows.Current_uA==round(STIM.sorted_current_uA(k)),1);
                if ~isempty(idx) && strcmp(rows.Keep{idx},'N'), ex(k)=true; end
            end
            valid = find(~ex);
            
            d = load(fullfile(T_PRM.STORAGE_DIR,rat,sprintf('%s_DATA.mat',rat)),'LF_Shank');
            mags = zeros(1,T_PRM.shank_no);
            for sh = 1:T_PRM.shank_no
                Vsh = squeeze(d.LF_Shank(sh).V{sess}(valid,:,xs:xe));
                flat = reshape(Vsh,[],size(Vsh,3));
                mags(sh) = abs(min(mean(flat,1)));
            end
            [~,bs] = max(mags);
            
            Vb = squeeze(d.LF_Shank(bs).V{sess}(valid,1:maxChannels,xs:xe));
            meanCh94 = squeeze(mean(Vb,1));
            VsumY = VsumY + meanCh94;
        end

        % Old rats
        for r = oldIdx
            rat = T_PRM.RATS{r};
            rows = tbl(strcmp(string(tbl.Rat), rat) & strcmp(tbl.Session,STIM.epochs{sess}), :);
            ex = false(1, numel(STIM.sorted_current_uA));
            for k = 1:numel(STIM.sorted_current_uA)
                idx = find(rows.Current_uA==round(STIM.sorted_current_uA(k)),1);
                if ~isempty(idx) && strcmp(rows.Keep{idx},'N'), ex(k)=true; end
            end
            valid = find(~ex);
            
            d = load(fullfile(T_PRM.STORAGE_DIR,rat,sprintf('%s_DATA.mat',rat)),'LF_Shank');
            mags = zeros(1,T_PRM.shank_no);
            for sh = 1:T_PRM.shank_no
                Vsh = squeeze(d.LF_Shank(sh).V{sess}(valid,:,xs:xe));
                flat = reshape(Vsh,[],size(Vsh,3));
                mags(sh) = abs(min(mean(flat,1)));
            end
            [~,bs] = max(mags);
            
            Vb = squeeze(d.LF_Shank(bs).V{sess}(valid,1:maxChannels,xs:xe));
            meanCh94 = squeeze(mean(Vb,1));
            VsumO = VsumO + meanCh94;
        end

        Vmat = Vsum / nRats;
        VmatY = VsumY / numel(youngIdx);
        VmatO = VsumO / numel(oldIdx);

        contrast_Vmat = real(Vmat .^ gamma);
        CmatY = real(VmatY.^gamma);
        CmatO = real(VmatO.^gamma);

        % All rats
        figure;
        imagesc(t, depths, contrast_Vmat);
        axis xy;
        set(gca,'YDir','reverse');
        colormap('jet'); colorbar;
        clim([-0.8 0.8]);
        hold on;
        line([0 0], [min(depths) max(depths)], 'Color','w','LineWidth',3);
        xlabel('Time (ms)');
        ylabel('Depth (\mum)');
        title(sprintf('Average Response to %s Stimulation', STIM.epochs{sess}));
        pubify_figure_axis_robust(14,14);
        img_name = sprintf('LFP_ALL_%s', STIM.epochs{sess});
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
        close;

        % Young vs Old side by side
        figure('Units','normalized','Position',[.1 .1 .8 .7]);
        tL = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        title(tL, sprintf('Average LF Response to %s Stimulation', STIM.epochs{sess}),...
            'Fontsize',16,'Fontweight','bold');

        ax1 = nexttile;
        imagesc(ax1, t, depths, CmatY);
        axis(ax1,'xy'); set(ax1,'YDir','reverse');
        colormap(ax1,'jet'); colorbar(ax1);
        clim(ax1,[-0.75 0.75]);
        line(ax1,[0 0],[min(depths) max(depths)],'Color','w','LineWidth',3);
        xlabel(ax1,'Time (ms)'); ylabel(ax1,'Depth (\mum)');
        title(ax1,'Young Rats');
        pubify_figure_axis_robust(14,14);

        ax2 = nexttile;
        imagesc(ax2, t, depths, CmatO);
        axis(ax2,'xy'); set(ax2,'YDir','reverse');
        colormap(ax2,'jet'); colorbar(ax2);
        clim(ax2,[-0.75 0.75]);
        line(ax2,[0 0],[min(depths) max(depths)],'Color','w','LineWidth',3);
        xlabel(ax2,'Time (ms)'); ylabel(ax2,'Depth (\mum)');
        title(ax2,'Old Rats');
        pubify_figure_axis_robust(14,14);
        
        img_name = sprintf('LFP_age_sep_%s', STIM.epochs{sess});
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
        close;
    end
end

end
