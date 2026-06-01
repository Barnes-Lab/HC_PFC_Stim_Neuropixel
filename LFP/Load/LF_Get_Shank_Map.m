function LF_Get_Shank_Map(T_PRM)
% Function to analyze PL and IL responses across shanks and create shank maps

%% Initial parameters
allStimTbl = readtable(fullfile(T_PRM.STORAGE_DIR,'ALL_STIM_TIMES_KEEP_CHECK.csv'),...
    'VariableNamingRule','preserve');

% Storage dir
ShankDir = fullfile(T_PRM.STORAGE_DIR,'ShankInfo');
if ~exist(ShankDir,"dir")
    mkdir(ShankDir);
end

% Get time parameters (0ms to 120ms window)
first = T_PRM.RATS{1};
S = load(fullfile(T_PRM.STORAGE_DIR, first, sprintf('%s_DATA.mat',first)),...
    'PRM','STIM','LF_Shank');
STIM = S.STIM;
ts = -S.PRM.time_extract_preStim*1e3;
te = S.PRM.time_extract_postStim*1e3;
t_all = linspace(ts,te,size(S.LF_Shank(1).V{1},3));
[~,t0_idx] = min(abs(t_all-0));
[~,t120_idx] = min(abs(t_all-120));
t = t_all(t0_idx:t120_idx);

% Analysis parameters
PL_MAX_DEPTH = 3550;
IL_MIN_DEPTH = 3750;
regions = {'PL', 'IL'};
nRats = length(T_PRM.RATS);
PL_color = [0.15, 0.5, 0.6];
IL_color = [0.88, 0.35, 0.11];

%% Getting data
R = struct(); detailedData = [];

for iRat = 1:nRats
    rat = T_PRM.RATS{iRat};
    fprintf('Processing Rat %s...\n', rat);
    ratField = sprintf('Rat%s', rat);
    
    % Load rat data once
    S_rat = load(fullfile(T_PRM.STORAGE_DIR,rat,sprintf('%s_DATA.mat',rat)), 'LF_Shank');
    
    for sess = 1:2
        sessField = sprintf('sess%d', sess);
        
        % Get valid stimuli
        rows = allStimTbl(strcmp(string(allStimTbl.Rat),rat) & ...
            strcmp(allStimTbl.Session,STIM.epochs{sess}),:);
        ex = false(1,numel(STIM.sorted_current_uA));
        for k=1:numel(STIM.sorted_current_uA)
            idx = find(rows.Current_uA==round(STIM.sorted_current_uA(k)),1);
            if ~isempty(idx) && strcmp(rows.keep{idx},'N'), ex(k)=true; end
        end
        valid = find(~ex);
        if isempty(valid), continue; end
        
        for shIdx = 1:4
            shankField = sprintf('Shank%d', shIdx);
            allDepths = S_rat.LF_Shank(shIdx).SortedDepths;
            Vb_all = S_rat.LF_Shank(shIdx).V{sess}(valid,:,t0_idx:t120_idx);
            nStim = size(Vb_all,1);
            
            for regIdx = 1:2
                region = regions{regIdx};
                channel_mask = (strcmp(region,'PL') & allDepths <= PL_MAX_DEPTH) | ...
                    (strcmp(region,'IL') & allDepths >= IL_MIN_DEPTH);
                
                if ~any(channel_mask), continue; end
                
                Vb = Vb_all(:, channel_mask, :);
                region_depths = allDepths(channel_mask);
                nChan = size(Vb,2);
                
                minVals = nan(nStim,1); minChans = nan(nStim,1);
                
                for iStim = 1:nStim
                    V_window = squeeze(Vb(iStim,:,:)); % Already 0-120ms window
                    [chan_minima, ~] = min(V_window, [], 2);
                    [global_min, global_chan_idx] = min(chan_minima);
                    
                    minVals(iStim) = global_min;
                    minChans(iStim) = global_chan_idx;
                    
                    % Store detailed data
                    detailedData = [detailedData; {rat, sess, shIdx, region, iStim, ...
                        global_min, global_chan_idx, region_depths(global_chan_idx)}];
                end
                
                % Store in structure
                R.(ratField).(sessField).(shankField).(region).minVals = minVals;
                R.(ratField).(sessField).(shankField).(region).minChans = minChans;
                R.(ratField).(sessField).(shankField).(region).depths = region_depths;
                R.(ratField).(sessField).(shankField).(region).traces = Vb;
            end
        end
    end
end

%% Saving CSV and calculating shanks
% Create detailed CSV
detailedTable = table(detailedData(:,1), detailedData(:,2), detailedData(:,3), ...
    detailedData(:,4), detailedData(:,5), detailedData(:,6), detailedData(:,7), ...
    detailedData(:,8), 'VariableNames', {'Rat', 'Session', 'Shank', 'Region', ...
    'Stim', 'MinValue', 'MinChannel', 'ChannelDepth'});
writetable(detailedTable, fullfile(T_PRM.STORAGE_DIR, 'DetailedStimData.csv'));

% Calculate averages
summaryData = [];
for iRat = 1:nRats
    rat = T_PRM.RATS{iRat};
    ratField = sprintf('Rat%s', rat);
    if ~isfield(R, ratField), continue; end
    
    for shIdx = 1:4
        PL_all = []; IL_all = [];
        for sess = 1:2
            sessField = sprintf('sess%d', sess);
            shankField = sprintf('Shank%d', shIdx);
            if isfield(R.(ratField), sessField) && isfield(R.(ratField).(sessField), shankField)
                if isfield(R.(ratField).(sessField).(shankField), 'PL')
                    PL_all = [PL_all; R.(ratField).(sessField).(shankField).PL.minVals];
                end
                if isfield(R.(ratField).(sessField).(shankField), 'IL')
                    IL_all = [IL_all; R.(ratField).(sessField).(shankField).IL.minVals];
                end
            end
        end
        summaryData = [summaryData; {rat, shIdx, mean(PL_all,'omitnan'), mean(IL_all,'omitnan')}];
    end
end

% Save summary and create shank map
summaryTable = table(summaryData(:,1), summaryData(:,2), summaryData(:,3), summaryData(:,4), ...
    'VariableNames', {'Rat', 'Shank', 'PL_AvgMin', 'IL_AvgMin'});
writetable(summaryTable, fullfile(ShankDir, 'ShankSummary.csv'));

% Select Shank for L2 adn L5 (will verify with histology)
shankMap = [];
for iRat = 1:nRats
    rat = T_PRM.RATS{iRat};
    ratData = summaryTable(strcmp(summaryTable.Rat, rat), :);
    if height(ratData) < 4, continue; end
    
    % L2 selection (shanks 1-2)
    s12 = ratData(ismember(cell2mat(ratData.Shank), [1,2]), :);
    [~, L2_idx] = min(cell2mat(s12.IL_AvgMin));
    if length(unique(cell2mat(s12.IL_AvgMin))) > 1 && abs(diff(cell2mat(s12.IL_AvgMin))) <= 0.25
        [~, L2_idx] = min(cell2mat(s12.PL_AvgMin));
    end
    
    % L5 selection (shanks 3-4)
    s34 = ratData(ismember(cell2mat(ratData.Shank), [3,4]), :);
    [~, L5_idx] = min(cell2mat(s34.IL_AvgMin));
    if length(unique(cell2mat(s34.IL_AvgMin))) > 1 && abs(diff(cell2mat(s34.IL_AvgMin))) <= 0.25
        [~, L5_idx] = min(cell2mat(s34.PL_AvgMin));
    end
    
    shankMap = [shankMap; {rat, cell2mat(s12.Shank(L2_idx)), cell2mat(s34.Shank(L5_idx))}];
end

writetable(table(shankMap(:,1), shankMap(:,2), shankMap(:,3), ...
    'VariableNames', {'Rat', 'L2_Shank', 'L5_Shank'}), ...
    fullfile(ShankDir, 'ShankMap.csv'));

%% Plotting
for iRat = 1:nRats
    rat = T_PRM.RATS{iRat};
    ratField = sprintf('Rat%s', rat);
    if ~isfield(R, ratField), continue; end
    
    % Initialize figure
    fig = figure('Position', [100, 100, 1600, 800]);
    sgtitle(sprintf('Rat %s', rat), 'FontSize', 16);
    colors = {PL_color, IL_color};
    
    for sess = 1:2
        sessField = sprintf('sess%d', sess);
        if ~isfield(R.(ratField), sessField), continue; end
        
        for shIdx = 1:4
            shankField = sprintf('Shank%d', shIdx);
            
            % Calculate subplot position: sess determines row, shIdx determines column
            subplot_pos = (sess-1)*4 + shIdx;
            subplot(2, 4, subplot_pos); 
            hold on;
            
            if ~isfield(R.(ratField).(sessField), shankField)
                title(sprintf('Sess %d - Shank %d (No Data)', sess, shIdx));
                continue;
            end
            
            for regIdx = 1:2
                region = regions{regIdx};
                if ~isfield(R.(ratField).(sessField).(shankField), region), continue; end
                
                data = R.(ratField).(sessField).(shankField).(region);
                
                % Get all traces and average them
                all_traces = [];
                for iStim = 1:length(data.minVals)
                    trace = squeeze(data.traces(iStim, data.minChans(iStim), :));
                    all_traces = [all_traces; trace'];
                end
                
                % Calculate average trace across all stimuli
                avg_trace = mean(all_traces, 1);
                
                color = colors{regIdx};
                plot(t, avg_trace, 'Color', color, 'LineWidth', 2);
            end
            
           
            title(sprintf('Sess %d - Shank %d', sess, shIdx));
            if sess == 2, xlabel('Time (ms)'); end
            if shIdx == 1, ylabel('Voltage'); end
            xlim([t(1), t(end)]);
        end
    end
    
    saveas(fig, fullfile(ShankDir, sprintf('MinTraces_Rat%s_BothSessions.png', rat)));
    close(fig);
end


fprintf('Analysis complete. Files saved in: %s\n', ShankDir);
end
