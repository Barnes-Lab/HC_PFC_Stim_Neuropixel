function NP2_Plot_LFP_DepthTraces(T_PRM)
% NP2_PLOT_LFP_DEPTHTRACES
%   Averages LFP over valid stims, remaps to a common 94‐bin depth axis,
%   then for each rat/session makes:
%     • Offset‐trace figure (1×2: L2/L5) with depth‐colored, inverted‐depth ordering
%     • Overlay‐trace figure (1×2) linking y‐axes and marking t=0
%   And for each session:
%     • Group‐avg offset figure (2×2: Y/O × L2/L5)
%     • Group‐avg overlay figure (2×2) linking y‐axes and marking t=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rat_plt=true;
grp_plt=true;

%Init params and dirs
smoothWin   = 1;
PRM_IMG_DIR  = fullfile(T_PRM.IMG_DIR,'DepthTraces');
if ~exist(PRM_IMG_DIR,'dir'), mkdir(PRM_IMG_DIR); end

overlay_dir=fullfile(PRM_IMG_DIR,'overlay_depths');
if ~exist(overlay_dir,'dir'), mkdir(overlay_dir); end
offset_dir=fullfile(PRM_IMG_DIR,'offset_depths');
if ~exist(offset_dir,'dir'), mkdir(offset_dir); end

% Load tables
allStimTbl = readtable(fullfile(T_PRM.STORAGE_DIR,'ALL_STIM_TIMES_KEEP_CHECK.csv'),...
    'VariableNamingRule','preserve');
shankMap   = readtable(fullfile(T_PRM.STORAGE_DIR,'ShankInfo','ShankMap.csv'),...
    'VariableNamingRule','preserve');

% Get times
first  = T_PRM.RATS{1};
S      = load(fullfile(T_PRM.STORAGE_DIR, first, sprintf('%s_DATA.mat',first)),...
    'PRM','STIM','LF_Shank');
PRM    = S.PRM; STIM = S.STIM;
ts     = -PRM.time_extract_preStim*1e3;
te     =  PRM.time_extract_postStim*1e3;
t_all  = linspace(ts,te,size(S.LF_Shank(1).V{1},3));
[~,xs] = min(abs(t_all+50)); [~,xe] = min(abs(t_all-150));
t      = t_all(xs:xe);
% Init params to loop over
shankTypes = {'L2','L5'};
ages       = {'Y','O'};
nRats      = numel(T_PRM.RATS);

%% Make a depth axis across shanks 94 chan
allDepths = [];
D0 = load(fullfile(T_PRM.STORAGE_DIR, first, sprintf('%s_DATA.mat',first)),'LF_Shank');
for sh=1:numel(D0.LF_Shank)
    allDepths = [allDepths; D0.LF_Shank(sh).SortedDepths];
end
globalDepths = linspace(min(allDepths), max(allDepths), 94)';

% depth‐color map (orange=shallow → purple=deep)
colors = [
    0.88,0.35,0.11;  % orange
    0.9,0.85,0.26;   % yellow
    0.55,0.7,0.4;    % light green
    0.15,0.5,0.6;    % teal
    0.25,0.41,0.88;  % blue
    0.35,0.05,0.7];  % purple
x   = linspace(0,1,size(colors,1));
xq  = linspace(0,1,94);
colorMap = interp1(x, colors, xq, 'linear');

% map each shank's sorted depths to global bins
mapRaw = cell(numel(D0.LF_Shank),1);
for sh=1:numel(D0.LF_Shank)
    dpths = D0.LF_Shank(sh).SortedDepths;
    mapRaw{sh} = arrayfun(@(d) find(abs(globalDepths-d)==min(abs(globalDepths-d)),1), dpths);
end
%% Loop over sessions
for sess = 1:2
    %initialize group collectors
    clear grpRaw
    for sType=1:2
        grpRaw.(shankTypes{sType}).Y = [];
        grpRaw.(shankTypes{sType}).O = [];
    end

    %% Per‐rat Plots
    for iRat = 1:nRats
        rat = T_PRM.RATS{iRat}; age = T_PRM.AGES{iRat};

        % valid stims
        rows = allStimTbl(strcmp(string(allStimTbl.Rat),rat)&...
            strcmp(allStimTbl.Session,STIM.epochs{sess}),:);
        ex = false(1,numel(STIM.sorted_current_uA));
        for k=1:numel(STIM.sorted_current_uA)
            idx = find(rows.Current_uA==round(STIM.sorted_current_uA(k)),1);
            if ~isempty(idx)&&strcmp(rows.keep{idx},'N'), ex(k)=true; end
        end
        valid = find(~ex);
        if isempty(valid), continue; end

        
            %% OFFSET‐TRACE FIGURE
            h1 = figure('Units','normalized','Position',[.1 .1 .9 .8],'Visible','on');
            tL1= tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
            title(tL1,sprintf('%s — %s mPFC DepthTraces (offset)',rat,STIM.epochs{sess}),...
                'FontSize',16,'FontWeight','bold');

            for sType=1:2
                shLabel = shankTypes{sType};
                rm    = strcmp(string(shankMap.Rat),string(rat));
                shIdx = shankMap.(shLabel)(rm);
                Vb    = squeeze(load(fullfile(T_PRM.STORAGE_DIR,rat,...
                    sprintf('%s_DATA.mat',rat)),'LF_Shank').LF_Shank(shIdx).V{sess}(valid,:,xs:xe));
                avgRaw= squeeze(mean(Vb,1));
                sRaw  = movmean(avgRaw,smoothWin,2);

                gRaw  = nan(94,numel(t));
                for ch=1:size(sRaw,1)
                    bin = mapRaw{shIdx}(ch);
                    gRaw(bin,:) = mean([gRaw(bin,:); sRaw(ch,:)],1,'omitnan');
                end
                grpRaw.(shLabel).(age)(end+1,:,:) = gRaw;

                ax = nexttile(sType); hold(ax,'on');
                offsetSpacing = 0.05;
                for bin=1:94
                    y = gRaw(bin,:) + (94-bin)*offsetSpacing;   % invert depth order
                    plot(ax,t,y,'LineWidth',1.5,'Color',colorMap(bin,:));
                end
                ax.YTick = [];  % remove y‐ticks
                xlabel(ax,'Time (ms)','FontSize',14);
                title(ax,sprintf('Layer %s',shLabel),'FontSize',14);
                xlim(ax,[min(t),max(t)]);
                line([0 0],ylim(ax),'Color','#A2142F','LineWidth',3);
                if sType==2
                    % color bar outside second plot
                    cbAx = axes('Position',[0.95 0.1 0.01 0.80],'YDir','reverse');
                    imagesc(cbAx, globalDepths);
                    colormap(cbAx, (colorMap));
                    yt = round(linspace(globalDepths(1),globalDepths(end),numel(cbAx.YTick))/100)*100;
                    set(cbAx,'YTickLabel',yt,'FontSize',10,'YAxisLocation','right','YTickLabelRotation',30);
                    ylabel(cbAx,'Depth (um)','FontSize',14);
                end
            end
            pubify_figure_axis_robust(14,14)
            saveas(h1,fullfile(offset_dir,...
                sprintf('%s_%s_DepthTraces_offset.png',rat,STIM.epochs{sess})));
            saveas(h1,fullfile(offset_dir,...
                sprintf('%s_%s_DepthTraces_offset.fig',rat,STIM.epochs{sess})));
            close(h1);

            %% OVERLAY‐TRACE FIGURE
            h2 = figure('Units','normalized','Position',[.1 .1 .8 .6],'Visible','off');
            tL2= tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
            title(tL2,sprintf('%s — %s Stim',rat,STIM.epochs{sess}),...
                'FontSize',16,'FontWeight','bold');

            axOL = gobjects(2,1);
            for sType=1:2
                shLabel = shankTypes{sType};
                axOL(sType) = nexttile(sType); hold(axOL(sType),'on');
                data = squeeze(grpRaw.(shLabel).(age)(end,:,:));  % this rat only
                for bin=1:94
                    plot(axOL(sType),t,data(bin,:),'LineWidth',1.5,'Color',colorMap(bin,:));
                end
                %axOL(sType).YTick = [];
                xlabel(axOL(sType),'Time (ms)','FontSize',14);
                ylabel(axOL(sType),'Voltage (mV)','FontSize',14);
                title(axOL(sType),sprintf('Layer %s',shLabel),'FontSize',14);
                xlim(axOL(sType),[-30,150]);
            end
            % link y‐axes and draw t=0 line
            linkaxes(axOL,'y');
            for k=1:2
                line([0 0], ylim(axOL(k)), 'Color','#A2142F','LineWidth',4, 'Parent', axOL(k));
            end
            % color bar
            cbAx2 = axes('Position',[0.95 0.1 0.01 0.80],'YDir','normal');
            imagesc(cbAx2, globalDepths);
            colormap(cbAx2, flip(colorMap));
            yt = round(linspace(globalDepths(1),globalDepths(end),numel(cbAx2.YTick))/100)*100;
            set(cbAx2,'YTickLabel',yt,'FontSize',10,'YAxisLocation','right','YTickLabelRotation',30);
            ylabel(cbAx2,'Depth (um)','FontSize',14);
            pubify_figure_axis_robust(14,14);
            saveas(h2,fullfile(overlay_dir,...
                sprintf('%s_%s_DepthTraces_overlay.png',rat,STIM.epochs{sess})));
            close(h2);
        

    end % per‐rat

    %% GROUP PLOTS
    if grp_plt
        %% — GROUP OFFSET: L2 (Young vs Old) —
        hOffL2 = figure('Units','normalized','Position',[.1 .1 .8 .6],'Visible','on');
        tOffL2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        title(tOffL2, sprintf('GroupAvg — %s mPFC L2 offset', STIM.epochs{sess}), ...
            'FontSize',16,'FontWeight','bold');

        axOffL2 = gobjects(2,1);
        offsetSpacing = 0.05;
        for a = 1:2
            axOffL2(a) = nexttile(a); hold(axOffL2(a),'on');
            data = squeeze(mean(grpRaw.L2.(ages{a}),1));      % [94×nT]
            sData = movmean(data, smoothWin, 2);
            for bin = 1:94
                y = sData(bin,:) + (94-bin)*offsetSpacing;   % invert depth
                plot(axOffL2(a), t, y, 'LineWidth',1.5, 'Color', colorMap(bin,:));
            end
            axOffL2(a).YTick = [];
            xlabel(axOffL2(a),'Time (ms)','FontSize',14);
            ylabel(axOffL2(a),'Voltage + offset','FontSize',14);
            title(axOffL2(a), ages{a}, 'FontSize',14);
            xlim(axOffL2(a), [-30,150]);
        end
        linkaxes(axOffL2,'y');
        for k = 1:2
            line([0 0], ylim(axOffL2(k)), 'Color','#A2142F','LineWidth',4, 'Parent', axOffL2(k));
        end
        % colorbar to the right of the second subplot
        cbAxOff2 = axes('Position',[0.95 0.1 0.01 0.80],'YDir','normal');
        imagesc(cbAxOff2, globalDepths);
        colormap(cbAxOff2, flip(colorMap));
        yt = round(linspace(globalDepths(1),globalDepths(end),numel(cbAxOff2.YTick))/100)*100;
        set(cbAxOff2,'YTickLabel',yt,'FontSize',10,'YAxisLocation','right','YTickLabelRotation',30);
        ylabel(cbAxOff2,'Depth (µm)','FontSize',14);
        pubify_figure_axis_robust(14,14);
        saveas(hOffL2, fullfile(offset_dir, sprintf('GroupAvg_%s_L2_offset.png', STIM.epochs{sess})));
        close(hOffL2);

        %% — GROUP OFFSET: L5 (Young vs Old) —
        hOffL5 = figure('Units','normalized','Position',[.1 .1 .8 .6],'Visible','on');
        tOffL5 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        title(tOffL5, sprintf('GroupAvg — %s mPFC L5 offset', STIM.epochs{sess}), ...
            'FontSize',16,'FontWeight','bold');

        axOffL5 = gobjects(2,1);
        for a = 1:2
            axOffL5(a) = nexttile(a); hold(axOffL5(a),'on');
            data = squeeze(mean(grpRaw.L5.(ages{a}),1));      % [94×nT]
            sData = movmean(data, smoothWin, 2);
            for bin = 1:94
                y = sData(bin,:) + (94-bin)*offsetSpacing;
                plot(axOffL5(a), t, y, 'LineWidth',1.5, 'Color', colorMap(bin,:));
            end
            axOffL5(a).YTick = [];
            xlabel(axOffL5(a),'Time (ms)','FontSize',14);
            ylabel(axOffL5(a),'Voltage + offset','FontSize',14);
            title(axOffL5(a), ages{a}, 'FontSize',14);
            xlim(axOffL5(a), [-30,150]);
        end
        linkaxes(axOffL5,'y');
        for k = 1:2
            line([0 0], ylim(axOffL5(k)), 'Color','#A2142F','LineWidth',4, 'Parent', axOffL5(k));
        end
        cbAxOff5 = axes('Position',[0.95 0.1 0.01 0.80],'YDir','normal');
        imagesc(cbAxOff5, globalDepths);
        colormap(cbAxOff5, flip(colorMap));
        yt = round(linspace(globalDepths(1),globalDepths(end),numel(cbAxOff5.YTick))/100)*100;
        set(cbAxOff5,'YTickLabel',yt,'FontSize',10,'YAxisLocation','right','YTickLabelRotation',30);
        ylabel(cbAxOff5,'Depth (µm)','FontSize',14);
        pubify_figure_axis_robust(14,14);
        saveas(hOffL5, fullfile(offset_dir, sprintf('GroupAvg_%s_L5_offset.png', STIM.epochs{sess})));
        close(hOffL5);
        %% GROUP OVERLAY: L2 (Young vs Old) —
        hL2 = figure('Units','normalized','Position',[.1 .1 .8 .6],'Visible','on');
        tL2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        title(tL2, sprintf('GroupAvg — %s mPFC L2 overlay', STIM.epochs{sess}), ...
            'FontSize',16,'FontWeight','bold');

        axL2 = gobjects(2,1);
        for a = 1:2
            axL2(a) = nexttile(a); hold(axL2(a),'on');
            data = squeeze(mean(grpRaw.L2.(ages{a}),1));    % [94×nT]
            for bin = 1:94
                plot(axL2(a), t, data(bin,:), ...
                    'LineWidth',1.5, 'Color', colorMap(bin,:));
            end
            if a==1
               ylabel(axL2(a),'Voltage (mV)','FontSize',14);
            end
            xlabel(axL2(a),'Time (ms)','FontSize',14);
            title(axL2(a), ages{a}, 'FontSize',14);
            xlim(axL2(a), [-30,120]);
        end
        % link y-axes and draw t=0
        linkaxes(axL2,'y');
        for k = 1:2
            line([0 0], ylim(axL2(k)), 'Color','#A2142F','LineWidth',4, 'Parent', axL2(k));
        end
        pubify_figure_axis_robust(14,14);
        saveas(hL2, fullfile(overlay_dir, sprintf('GroupAvg_%s_L2_overlay_new.png', STIM.epochs{sess})));
        close(hL2);

        %% — GROUP OVERLAY: L5 (Young vs Old) —
        hL5 = figure('Units','normalized','Position',[.1 .1 .8 .6],'Visible','on');
        tL5 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        title(tL5, sprintf('GroupAvg — %s mPFC L5 overlay', STIM.epochs{sess}), ...
            'FontSize',16,'FontWeight','bold');

        axL5 = gobjects(2,1);
        for a = 1:2
            axL5(a) = nexttile(a); hold(axL5(a),'on');
            data = squeeze(mean(grpRaw.L5.(ages{a}),1));    % [94×nT]
            for bin = 1:94
                plot(axL5(a), t, data(bin,:), ...
                    'LineWidth',1.5, 'Color', colorMap(bin,:));
            end
            if a==1
               ylabel(axL5(a),'Voltage (mV)','FontSize',14);
            end
            xlabel(axL5(a),'Time (ms)','FontSize',14);
            title(axL5(a), ages{a}, 'FontSize',14);
            xlim(axL5(a), [-30,120]);
        end
        % link y-axes and draw t=0
        linkaxes(axL5,'y');
        for k = 1:2
            line([0 0], ylim(axL5(k)), 'Color','#A2142F','LineWidth',3, 'Parent', axL5(k));
        end
           pubify_figure_axis_robust(14,14);
        saveas(hL5, fullfile(overlay_dir, sprintf('GroupAvg_%s_L5_overlay_new.png', STIM.epochs{sess})));
        close(hL5);
    end

    % Calculate session means 
    grpRawAll{sess}.L2.Y = grpRaw.L2.Y;   % size = [nRats × 94 × nT]
    grpRawAll{sess}.L2.O = grpRaw.L2.O;
    grpRawAll{sess}.L5.Y = grpRaw.L5.Y;
    grpRawAll{sess}.L5.O = grpRaw.L5.O;
end
% Combine iHC and vHC for each group/layer
groupMean = struct();

for layer = {'L2','L5'}
    L = layer{1};
    for a = 1:2
        age = ages{a}; % 'Y' or 'O'
        data1 = squeeze(mean(grpRawAll{1}.(L).(age), 1)); % [nCh x nT], iHC
        data2 = squeeze(mean(grpRawAll{2}.(L).(age), 1)); % [nCh x nT], vHC

        groupMean.(L).(age) = (data1 + data2) / 2;
    end
end
% Plot overall offsets
for li = 1:2
    layer = {'L2','L5'};
    L = layer{li};
    figure('Units','normalized','Position',[.1 .1 .8 .6]);
    tL = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    title(tL, sprintf('Average mPFC LFP: iHC + vHC %s',L), 'FontSize', 16, 'FontWeight', 'bold');
    for a = 1:2
        ax(a) = nexttile(a); hold(ax(a),'on');
        data = groupMean.(L).(ages{a});

        for bin = 1:94
            plot(ax(a), t, data(bin,:), ...
                'LineWidth',1.5, 'Color', colorMap(bin,:));
        end
        if a==1
            ylabel(ax(a),'Voltage (mV)','FontSize',14);
        end
        xlabel(ax(a),'Time (ms)','FontSize',14);
        title(ax(a), ages{a}, 'FontSize',14);
        xlim(ax(a), [-20,120]);
        pubify_figure_axis_robust(16,16);
    end
    % link y-axes and draw t=0
    linkaxes(ax,'y');
    for k = 1:2
        line([0 0], ylim(ax(k)), 'Color','#A2142F','LineWidth',3, 'Parent', ax(k));
    end
    
    saveas(gcf, fullfile(overlay_dir, sprintf('AllSess_%s.png',L)));
    %close(gcf);
end



end
