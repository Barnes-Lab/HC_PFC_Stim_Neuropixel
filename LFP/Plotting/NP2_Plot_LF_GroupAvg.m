function NP2_Plot_LF_GroupAvg(T_PRM)
% NP2_PLOT_LF_GROUPAVG
%   Averages IL/PL (mapped shank) across valid stims, rats & sessions.
%   Creates per-rat traces (L2 left, L5 right with IL/PL overlaid) and combined session averages.

%% 0) Load metadata & time vector
PRM_IMG_DIR = fullfile(T_PRM.IMG_DIR,"RegionAverages");
if ~exist(PRM_IMG_DIR,"dir"), mkdir(PRM_IMG_DIR); end

allStimTbl = readtable(fullfile(T_PRM.STORAGE_DIR,'ALL_STIM_TIMES_KEEP_CHECK.csv'),...
    'VariableNamingRule','preserve');
first = T_PRM.RATS{1};
load(fullfile(T_PRM.STORAGE_DIR, first, sprintf('%s_DATA.mat',first)),...
    'PRM','LF_Shank','STIM');

ts = -PRM.time_extract_preStim*1e3;
te =  PRM.time_extract_postStim*1e3;
t_all = linspace(ts,te,size(LF_Shank(1).V{1},3));
[~,xs] = min(abs(t_all+50));
[~,xe] = min(abs(t_all-150));
t = t_all(xs:xe);
PRM_IL_PL_gap = 4;

%% Load shank map & prepare
shankMap = readtable(fullfile(T_PRM.STORAGE_DIR,'ShankInfo','ShankMap.csv'),...
    'VariableNamingRule','preserve');
shankTypes = {'L2','L5'};

ages = {'Y','O'};
nRats = numel(T_PRM.RATS);

% CONSISTENT COLORS: IL=orange, PL=teal
baseColor.IL = [0.88, 0.35, 0.11];  % Orange
baseColor.PL = [0.15, 0.5, 0.6];     % Teal
darkFactor = 0.5;
darkColor.IL = baseColor.IL * darkFactor;
darkColor.PL = baseColor.PL * darkFactor;
lineStyles = {'-','--'};  % Y solid, O dashed

%% Extract data for all shanks first, then create per-rat figures
for sess = 1:2
    for iRat = 1:nRats
        rat = T_PRM.RATS{iRat};
        
        % Storage for this rat's data
        ratData = struct();
        
        % Loop over L2 and L5 to extract data
        for sType = 1:2
            shankLabel = shankTypes{sType};
            
            % Valid stims
            rows = allStimTbl(strcmp(string(allStimTbl.Rat),rat) & ...
                strcmp(allStimTbl.Session,STIM.epochs{sess}), :);
            ex = false(1, numel(STIM.sorted_current_uA));
            for k = 1:numel(STIM.sorted_current_uA)
                idx = find(rows.Current_uA==round(STIM.sorted_current_uA(k)),1);
                if ~isempty(idx) && strcmp(rows.keep{idx},'N')
                    ex(k) = true;
                end
            end
            valid = find(~ex);
            
            % Mapped shank lookup
            rowMap = strcmp(string(shankMap.Rat), rat);
            if ~any(rowMap)
                error('Rat %s missing from ShankMap.csv', rat);
            end
            shank_idx = shankMap.(shankLabel)(rowMap);
            
            % Load data & slice
            D = load(fullfile(T_PRM.STORAGE_DIR,rat,sprintf('%s_DATA.mat',rat)), 'LF_Shank');
            Vb = squeeze(D.LF_Shank(shank_idx).V{sess}(valid,:,xs:xe));
            
            % Define IL & PL channels
            ilc = 1 : D.LF_Shank(shank_idx).IL_end_chan;
            plc = ilc(end)+PRM_IL_PL_gap : size(Vb,2);
            
            % IL region filtering & mean
            dat = Vb(:, ilc, :);
            pm = zeros(size(dat,1), numel(t));
            for s = 1:size(dat,1)
                m = squeeze(dat(s,:,:));
                keep = min(m,[],2) <= min(m(:))*0.5;
                if ~any(keep), keep = true(size(keep)); end
                pm(s,:) = mean(m(keep,:),1);
            end
            ratData.(shankLabel).IL = mean(pm,1);
            
            % PL region filtering & mean
            dat = Vb(:, plc, :);
            pm = zeros(size(dat,1), numel(t));
            for s = 1:size(dat,1)
                m = squeeze(dat(s,:,:));
                keep = min(m,[],2) <= min(m(:))*0.5;
                if ~any(keep), keep = true(size(keep)); end
                pm(s,:) = mean(m(keep,:),1);
            end
            ratData.(shankLabel).PL = mean(pm,1);
        end
        
        %% PER-RAT FIGURE: L2 on left, L5 on right, IL and PL overlaid in each
        figure('Units','normalized','Position',[.1 .1 .8 .5]);
        tL = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        title(tL, sprintf('%s - %s', rat, STIM.epochs{sess}),...
            'FontSize',16,'FontWeight','bold');
        
        for sType = 1:2
            shankLabel = shankTypes{sType};
            ax = nexttile;
            hold(ax,'on');
            
            % Plot IL and PL overlaid
            hIL = plot(ax, t, ratData.(shankLabel).IL, 'LineWidth',2, 'Color',baseColor.IL);
            hPL = plot(ax, t, ratData.(shankLabel).PL, 'LineWidth',2, 'Color',baseColor.PL);
            line(ax, [0 0], ylim(ax), 'Color','#A2142F','LineWidth',2);
            
            xlabel(ax,'Time (ms)'); ylabel(ax,'Voltage (mV)');
            title(ax,shankLabel);
            xlim(ax,[-50 150]);
            legend(ax,[hIL hPL],{'IL','PL'},'Location','northeast');
            pubify_figure_axis_robust(14,14);
        end
        
        img_name = sprintf('%s_%s_L2_L5_ILPL', rat, STIM.epochs{sess});
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
        close;
    end
end

%% Now do group analysis
for sType = 1:2
    shankLabel = shankTypes{sType};
    fprintf('\n===== Running for shank %s =====\n', shankLabel);

    % Storage for combined sessions
    allIL_byAgeRat_combined = cell(2,1);
    allPL_byAgeRat_combined = cell(2,1);
    
    % GROUP-AVG extraction per session
    allIL_byAgeRat = cell(2,2);
    allPL_byAgeRat = cell(2,2);
    clear grp
    for a = 1:2
        rats_a = find(strcmp(T_PRM.AGES, ages{a}));
        nR = numel(rats_a);
        
        % Combined storage
        ratIL_combined = zeros(nR, numel(t));
        ratPL_combined = zeros(nR, numel(t));
        
        for sess = 1:2
            ratIL = zeros(nR, numel(t));
            ratPL = zeros(nR, numel(t));

            for j = 1:nR
                iRat = rats_a(j);
                rat = T_PRM.RATS{iRat};

                % Valid stims
                rows = allStimTbl(strcmp(string(allStimTbl.Rat),rat) & ...
                    strcmp(allStimTbl.Session,STIM.epochs{sess}), :);
                ex = false(1, numel(STIM.sorted_current_uA));
                for k = 1:numel(STIM.sorted_current_uA)
                    idx = find(rows.Current_uA==round(STIM.sorted_current_uA(k)),1);
                    if ~isempty(idx) && strcmp(rows.keep{idx},'N')
                        ex(k) = true;
                    end
                end
                valid = find(~ex);

                % Mapped shank lookup
                rowMap = strcmp(string(shankMap.Rat), rat);
                if ~any(rowMap)
                    error('Rat %s missing from ShankMap.csv', rat);
                end
                shank_idx = shankMap.(shankLabel)(rowMap);

                % Load data & slice
                D = load(fullfile(T_PRM.STORAGE_DIR,rat,sprintf('%s_DATA.mat',rat)), 'LF_Shank');
                Vb = squeeze(D.LF_Shank(shank_idx).V{sess}(valid,:,xs:xe));

                % Define IL & PL channels
                ilc = 1 : D.LF_Shank(shank_idx).IL_end_chan;
                plc = ilc(end)+PRM_IL_PL_gap : size(Vb,2);

                % IL region filtering & mean
                dat = Vb(:, ilc, :);
                pm = zeros(size(dat,1), numel(t));
                for s = 1:size(dat,1)
                    m = squeeze(dat(s,:,:));
                    keep = min(m,[],2) <= min(m(:))*0.5;
                    if ~any(keep), keep = true(size(keep)); end
                    pm(s,:) = mean(m(keep,:),1);
                end
                ratIL(j,:) = mean(pm,1);

                % PL region filtering & mean
                dat = Vb(:, plc, :);
                pm = zeros(size(dat,1), numel(t));
                for s = 1:size(dat,1)
                    m = squeeze(dat(s,:,:));
                    keep = min(m,[],2) <= min(m(:))*0.5;
                    if ~any(keep), keep = true(size(keep)); end
                    pm(s,:) = mean(m(keep,:),1);
                end
                ratPL(j,:) = mean(pm,1);
            end

            allIL_byAgeRat{a,sess} = ratIL;
            allPL_byAgeRat{a,sess} = ratPL;

            grp(a,sess).meanIL = mean(ratIL,1);
            grp(a,sess).semIL = std(ratIL,0,1)./sqrt(nR);
            grp(a,sess).meanPL = mean(ratPL,1);
            grp(a,sess).semPL = std(ratPL,0,1)./sqrt(nR);
            
            % Accumulate for combined average
            ratIL_combined = ratIL_combined + ratIL;
            ratPL_combined = ratPL_combined + ratPL;
        end
        
        % Average across sessions
        ratIL_combined = ratIL_combined / 2;
        ratPL_combined = ratPL_combined / 2;
        
        allIL_byAgeRat_combined{a} = ratIL_combined;
        allPL_byAgeRat_combined{a} = ratPL_combined;
    end

    %% Plot combined session average (iHC + vHC) for this shank
    figure; hold on;
    hLines = gobjects(4,1);
    labels = strings(4,1);
    cnt = 1;
    for a = 1:2
        ls = lineStyles{a};
        ratIL = allIL_byAgeRat_combined{a};
        ratPL = allPL_byAgeRat_combined{a};
        
        for R = ["IL","PL"]
            C = ternary(a==1, baseColor.(R), darkColor.(R));
            if R=="IL"
                mn = mean(ratIL,1);
                se = std(ratIL,0,1)./sqrt(size(ratIL,1));
            else
                mn = mean(ratPL,1);
                se = std(ratPL,0,1)./sqrt(size(ratPL,1));
            end
            patch([t, fliplr(t)], [mn+se, fliplr(mn-se)], C, ...
                'FaceAlpha',0.2,'EdgeColor','none');
            hLines(cnt) = plot(t, mn, 'LineStyle',ls, 'LineWidth',2,'Color',C);
            labels(cnt) = sprintf('%s %s (%s)', ages{a}, R, shankLabel);
            cnt = cnt+1;
        end
    end
    xlabel('Time (ms)'); ylabel('Voltage (mV)');
    title(sprintf('Combined Sessions (iHC+vHC) — %s', shankLabel));
    line([0 0], ylim, 'Color','#A2142F','LineWidth',2);
    legend(hLines, labels,'Location','northeast');
    xlim([-50 150]);
    pubify_figure_axis_robust(14,14);
    img_name = sprintf('GroupAvg_Combined_Sessions_%s', shankLabel);
    saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
    saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
    close;

    %% Plot group avg ±SEM per session
    for sess = 1:2
        figure; hold on;
        hLines = gobjects(4,1);
        labels = strings(4,1);
        cnt = 1;
        for a = 1:2
            ls = lineStyles{a};
            for R = ["IL","PL"]
                C = ternary(a==1, baseColor.(R), darkColor.(R));
                mn = grp(a,sess).("mean"+R);
                se = grp(a,sess).("sem"+R);
                patch([t, fliplr(t)], [mn+se, fliplr(mn-se)], C, ...
                    'FaceAlpha',0.2,'EdgeColor','none');
                hLines(cnt) = plot(t, mn, 'LineStyle',ls, 'LineWidth',2,'Color',C);
                labels(cnt) = sprintf('%s %s (%s)', ages{a}, R, shankLabel);
                cnt = cnt+1;
            end
        end
        xlabel('Time (ms)'); ylabel('Voltage (mV)');
        title(sprintf('Group avg ±SEM — %s — %s', STIM.epochs{sess}, shankLabel));
        line([0 0], ylim, 'Color','#A2142F','LineWidth',2);
        legend(hLines, labels,'Location','northeast');
        xlim([-50 150]);
        pubify_figure_axis_robust(14,14);
        img_name = sprintf('GroupAvg_%s_%s', STIM.epochs{sess}, shankLabel);
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
        close;
    end

    %% Per-rat traces + group mean overlay
    for sess = 1:2
        figure; hold on;
        hMean = gobjects(4,1); labs = strings(4,1); cnt = 1;
        for a = 1:2
            ls = lineStyles{a};
            for R = ["IL","PL"]
                C = ternary(a==1, baseColor.(R), darkColor.(R));
                dataMat = allIL_byAgeRat{a,sess};
                if R=="PL", dataMat = allPL_byAgeRat{a,sess}; end
                for iR = 1:size(dataMat,1)
                    h = plot(t, dataMat(iR,:), 'LineStyle',ls, 'LineWidth',2,'Color',C);
                    h.Color(4) = 0.2;
                end
                mn = grp(a,sess).("mean"+R);
                hMean(cnt) = plot(t, mn, 'LineStyle',ls, 'LineWidth',2,'Color',C);
                labs(cnt) = sprintf('%s %s (%s)', ages{a}, R, shankLabel);
                cnt = cnt+1;
            end
        end
        xlabel('Time (ms)'); ylabel('Voltage (mV)');
        title(sprintf('Per-Rat & Group Mean — %s — %s', STIM.epochs{sess}, shankLabel));
        line([0 0], ylim, 'Color','#A2142F','LineWidth',2);
        legend(hMean, labs,'Location','northeast');
        xlim([-50 150]);
        pubify_figure_axis_robust(14,14);
        img_name = sprintf('GroupAvg_AllRats_%s_%s', STIM.epochs{sess}, shankLabel);
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
        saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
        close;
    end
end

%% TOP-10 channel selection per age & session
topIL = cell(2,2);
topPL = cell(2,2);
for sess = 1:2
    for a = 1:2
        M_IL = []; M_PL = [];
        for iR = find(strcmp(T_PRM.AGES, ages{a}))
            rat = T_PRM.RATS{iR};
            rows = allStimTbl(strcmp(string(allStimTbl.Rat),rat) & ...
                strcmp(allStimTbl.Session,STIM.epochs{sess}),:);
            ex = false(1,numel(STIM.sorted_current_uA));
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
            [~, bs] = max(mags);
            
            Vb = squeeze(d.LF_Shank(bs).V{sess}(valid,:,xs:xe));
            ilc = 1:d.LF_Shank(bs).IL_end_chan;
            plc = ilc(end)+PRM_IL_PL_gap : size(Vb,2);
            M_IL = [M_IL; squeeze(mean(Vb(:,ilc,:),1))];
            M_PL = [M_PL; squeeze(mean(Vb(:,plc,:),1))];
        end
        
        minIL = min(M_IL, [], 2);
        [~, idxIL] = sort(minIL, 'ascend');
        nIL = min(10, numel(idxIL));
        topIL{a,sess} = M_IL(idxIL(1:nIL), :);

        minPL = min(M_PL, [], 2);
        [~, idxPL] = sort(minPL, 'ascend');
        nPL = min(10, numel(idxPL));
        topPL{a,sess} = M_PL(idxPL(1:nPL), :);
    end
end

%% PER-RAT mean abs-min across top-10 channels
minValsPL = nan(2,nRats);
minValsIL = nan(2,nRats);
for sess = 1:2
    for iR = 1:nRats
        rat = T_PRM.RATS{iR};
        rows = allStimTbl(strcmp(string(allStimTbl.Rat),rat) & ...
            strcmp(allStimTbl.Session,STIM.epochs{sess}),:);
        ex = false(1,numel(STIM.sorted_current_uA));
        for k = 1:numel(STIM.sorted_current_uA)
            idx = find(rows.Current_uA==round(STIM.sorted_current_uA(k)),1);
            if ~isempty(idx)&&strcmp(rows.Keep{idx},'N'), ex(k)=true; end
        end
        valid = find(~ex);
        d = load(fullfile(T_PRM.STORAGE_DIR,rat,sprintf('%s_DATA.mat',rat)),'LF_Shank');
        mags = zeros(1,T_PRM.shank_no);
        for sh = 1:T_PRM.shank_no
            Vsh = squeeze(d.LF_Shank(sh).V{sess}(valid,:,xs:xe));
            flat = reshape(Vsh,[],size(Vsh,3));
            mags(sh) = abs(min(mean(flat,1)));
        end
        [~, bs] = max(mags);
        Vb = squeeze(d.LF_Shank(bs).V{sess}(valid,:,xs:xe));
        ilc = 1:d.LF_Shank(bs).IL_end_chan;
        plc = ilc(end)+PRM_IL_PL_gap : size(Vb,2);
        
        chIL = squeeze(mean(Vb(:,ilc,:),1));
        mIL = abs(min(chIL,[],2));
        [~, ordIL] = sort(mIL,'descend');
        mIL10 = mIL(ordIL(1:min(10,numel(ordIL))));
        minValsIL(sess,iR) = mean(mIL10);
        
        chPL = squeeze(mean(Vb(:,plc,:),1));
        mPL = abs(min(chPL,[],2));
        [~, ordPL] = sort(mPL,'descend');
        mPL10 = mPL(ordPL(1:min(10,numel(ordPL))));
        minValsPL(sess,iR) = mean(mPL10);
    end
end

%% 2×2 SUMMARY FIGURE
figure('Units','normalized','Position',[.1 .1 .8 .7]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for sess = 1:2
    ax1 = nexttile(sess);
    hold(ax1,'on');
    hL = gobjects(4,1); labs = strings(4,1); cnt = 1;
    for a = 1:2
        for reg = {'IL','PL'}
            R = reg{1};
            C = ternary(a==1, baseColor.(R), darkColor.(R));
            ls = lineStyles{a};
            M = topIL{a,sess}; if strcmp(R,'PL'), M = topPL{a,sess}; end
            mu = mean(M,1); sd = std(M,0,1);
            patch(ax1, [t fliplr(t)], [mu+sd fliplr(mu-sd)], C, ...
                'FaceAlpha',.2, 'EdgeColor','none');
            hL(cnt) = plot(ax1, t, mu, 'LineStyle',ls, 'LineWidth',2, 'Color',C);
            labs(cnt) = sprintf('%s %s', ages{a}, R);
            cnt = cnt+1;
        end
    end
    xlabel(ax1,'Time (ms)'); ylabel(ax1,'Voltage (mV)');
    title(ax1, sprintf('Sess %d Top-10', sess));
    legend(ax1, hL, labs, 'Location','northeast');
    line(ax1, [0 0], ylim(ax1), 'Color','#A2142F','LineWidth',2);
    pubify_figure_axis_robust(14,14);

    ax2 = nexttile(sess+2);
    hold(ax2,'on');

    isY = strcmp(T_PRM.AGES,'Y');
    isO = strcmp(T_PRM.AGES,'O');
    Y_IL = mean(minValsIL(sess,isY));
    O_IL = mean(minValsIL(sess,isO));
    Y_PL = mean(minValsPL(sess,isY));
    O_PL = mean(minValsPL(sess,isO));
    barVals = [Y_IL, O_IL, Y_PL, O_PL];

    b = bar(ax2, 1:4, barVals, 'FaceColor','flat');
    b.CData = [baseColor.IL; darkColor.IL; baseColor.PL; darkColor.PL];

    for k = 1:nRats
        if isY(k)
            scatter(ax2, 1 + (rand*0.3 - 0.15), minValsIL(sess,k), 30, 'k','filled','MarkerEdgeColor','w');
            scatter(ax2, 3 + (rand*0.3 - 0.15), minValsPL(sess,k), 30, 'k','filled','MarkerEdgeColor','w');
        else
            scatter(ax2, 2 + (rand*0.3 - 0.15), minValsIL(sess,k), 30, 'k','filled','MarkerEdgeColor','w');
            scatter(ax2, 4 + (rand*0.3 - 0.15), minValsPL(sess,k), 30, 'k','filled','MarkerEdgeColor','w');
        end
    end

    ax2.XTick = 1:4;
    ax2.XTickLabel = {'Y-IL','O-IL','Y-PL','O-PL'};
    ylabel(ax2,'Mean abs-min');
    title(ax2, sprintf('Sess %d', sess));
    pubify_figure_axis_robust(14,14);
end

sgtitle('Absolute Minimum Response by Age & Region','FontSize',16);
img_name = 'Summary_Top10_and_Bars';
saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.png']));
saveas(gcf, fullfile(PRM_IMG_DIR, [img_name '.fig']));
close;
end
