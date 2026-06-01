function[]= NP2_Plot_AllStim(rat_no, data_filepath, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function[]= NP2_Plot_AllStim(rat_no, data_filepath, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%   rat_no:         String
%   data_filepath:  Path to the .mat file containing PRM, STIM, and LF_Shank
%   varargin:       Optional flags (same as NP2_Plot_LF) to enable/disable sections
%
% This function generates plots per session:
%  1) All concatenated stim traces per channel, stacked by shank (with exclusion of noisy stims)
%  2) Mean IL vs PL per shank, stacked by shank (with STD shading, exclusion of noisy stims, and red zero‐stim lines)
%  3) Mean IL vs PL averaged across shanks (with STD shading, exclusion of noisy stims, and red zero‐stim lines)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS TO CHANGE
PRM_IMG_DIR = 'C:\DATA\Neuropixel\NP2\Processed\LFP_RawTraces';
PRM_STORAGE_DIR='';
PL_line         = true;
PL_regions_line = true;
PL_region_avg   = true;
PRM_IL_PL_gap   = 2;              % Gap between IL and PL channel indices
IL_col = [0.35, 0.05, 0.7];       % Infralimbic color
PL_col = [0.88, 0.35, 0.11];      % Prelimbic color

% Plot window relative to stim (ms)
plot_preStim  = -20;
plot_postStim = 120;

Extract_varargin;  % override defaults if provided

%% Load data
load(data_filepath, 'PRM', 'STIM', 'LF_Shank');
stim_currents = STIM.sorted_current_uA;   % [1 x nStims]
nStims        = length(stim_currents);
nShanks       = PRM.shank_no;             % should be 4

% Time axis for one stim window (ms)
time_start = -PRM.time_extract_preStim * 1000;
time_end   = PRM.time_extract_postStim * 1000;
x_time_all = linspace(time_start, time_end, size(LF_Shank(1).V{1,1}, 3));  % full time axis
% Find indices corresponding to cropped window
[~, x_start] = min(abs(x_time_all - plot_preStim));
[~, x_end]   = min(abs(x_time_all - plot_postStim));
x_time = x_time_all(x_start:x_end);  % cropped time axis for plotting
nPts   = length(x_time);
window_duration = x_time(end) - x_time(1);  % in ms, duration of one stim window

% Precompute prestim indices and threshold
prestim_idx = find(x_time_all < 0);               % samples before 0 ms
dt_ms       = x_time_all(2) - x_time_all(1);       % sampling interval in ms
threshold_samples = round(5 / dt_ms);           % # samples corresponding to 5 seconds

% Base colormap for concatenated channels (for PL_line section)
base_colors = [
    0.35, 0.05, 0.7;   % Purple
    0.25, 0.41, 0.88;  % Blue
    0.15, 0.5, 0.6;    % Teal
    0.55, 0.7, 0.4;    % Light green
    0.9, 0.85, 0.26;   % Yellow
    0.88, 0.35, 0.11   % Orange
    ];

%% Identify and exclude noisy stim amplitudes per session
% === Read the CSV once at the top of the function ===
allStimTbl = readtable(fullfile(PRM_STORAGE_DIR, 'ALL_STIM_TIMES_KEEP_CHECK.csv'));

% For each session, build two valid‐stim lists: “check” and “keep” ===
validStimIdx_check = cell(1,2);
validStimIdx_keep  = cell(1,2);

for sess = 1:2
    session_name = STIM.epochs{sess};
    %  Extract only the rows for this rat & this session
    rows = allStimTbl( ...
        strcmp(string(allStimTbl.Rat), string(rat_no)) & ...
        strcmp(allStimTbl.Session, string(session_name)), : );
    currents_session = rows.Current_uA;  % numeric vector of currents in CSV
    
    excluded_check = false(1, nStims);
    excluded_keep  = false(1, nStims);
    
    for iStim = 1:nStims
        % Find the matching row in the CSV for this current
        idx = find(currents_session == round(stim_currents(iStim)), 1);
        if ~isempty(idx)
            % If Check=='Y', exclude from “check” set
            if strcmp(rows.Check{idx}, 'Y')
                excluded_check(iStim) = true;
            end
            % If Keep=='N', exclude from “keep” set
            if strcmp(rows.Keep{idx}, 'N')
                excluded_keep(iStim) = true;
            end
        end
    end
    
    validStimIdx_check{sess} = find(~excluded_check);
    validStimIdx_keep{sess}  = find(~excluded_keep);
end

%% All concatenated stim traces per channel, stacked by shank 
if PL_line
    for sess = 1:2
        session_name = STIM.epochs{sess};
        validStimIdx = validStimIdx_keep{sess};
        nValidStims  = numel(validStimIdx);
        if nValidStims == 0
            warning('Session %s: All stimuli excluded due to noise. Skipping concatenated plot.', session_name);
            continue;
        end

        % Compute y-axis limit across all shanks and valid stimuli
        minValSession = inf;
        maxValSession = -inf;
        for shank = 1:nShanks
            nChan = length(LF_Shank(shank).SortedDepths);
            for chan = 1:nChan
                for iVal = 1:nValidStims
                    iStim = validStimIdx(iVal);
                    dt = squeeze(LF_Shank(shank).V{sess}(iStim, chan, x_start:x_end));
                    maxValSession = max(maxValSession, max(dt(:)));
                    minValSession = min(minValSession, min(dt(:)));
                end
            end
        end
        % Then set y-limits using raw min/max (no abs)
        ylim([minValSession, maxValSession]);

        % Build concatenated x-axis offsets for valid stimuli
        stim_offsets = (0:(nValidStims-1)) * window_duration;  % in ms

        % Create figure with tiled layout: 4 rows × 1 column
        figure('Position', [50, 50, 2300, 1400],'Visible','off');
        tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        sgtitle(sprintf('AllStim %s - %s', rat_no, session_name), 'FontSize', 20);

        for shank = 1:nShanks
            nChan = length(LF_Shank(shank).SortedDepths);
            % Build interpolated colormap per channel
            xcol = linspace(0,1,size(base_colors,1));
            xq   = linspace(0,1,nChan);
            colorMap = zeros(nChan,3);
            for c = 1:3
                colorMap(:,c) = interp1(xcol, base_colors(:,c), xq, 'linear');
            end

            % Build concatenated x-axis
            x_concat = nan(1, nValidStims * nPts);
            for iVal = 1:nValidStims
                offset = stim_offsets(iVal);
                idxs   = (iVal-1)*nPts + (1:nPts);
                x_concat(idxs) = x_time + offset;
            end

            nexttile;
            hold on;
            ax = gca;
            title(sprintf('Shank %d', shank), 'FontSize', 16);
            if shank == nShanks
                xlabel('Stimulation Amplitude (\muA)', 'FontSize', 14);
                xticks(stim_offsets);
                xticklabels(arrayfun(@(c) sprintf('%.0f', stim_currents(c)), validStimIdx, 'UniformOutput', false));
            else
                ax.XTick = [];
            end
            ylabel('Voltage (mV)', 'FontSize', 14);
            
            set(gca, 'FontSize', 12);

            % Plot each channel's concatenated trace
            for chan = 1:nChan
                signal_concat = nan(1, nValidStims * nPts);
                for iVal = 1:nValidStims
                    iStim = validStimIdx(iVal);
                    dt = squeeze(LF_Shank(shank).V{sess}(iStim, chan, x_start:x_end));
                    signal_concat((iVal-1)*nPts + (1:nPts)) = dt;
                end
                plot(x_concat, signal_concat, 'Color', colorMap(chan,:), 'LineWidth', 0.5);
            end

            % Plot vertical stimulus lines and relative zero‐stim lines for each stim
            yl = ylim;
            for iVal = 1:nValidStims
                offset = stim_offsets(iVal);
                % Black line at start of each window:
                plot([offset offset], yl, 'Color', 'k', 'LineWidth', 1);
                plot([offset + window_duration offset + window_duration], yl, 'Color', 'k', 'LineWidth', 1);
                % Red line at relative 0 ms within each stim window:
                plot([offset offset], yl, 'Color', 'r', 'LineWidth', 1.5);
            end

            hold off;
        end
        cbAx = axes('Position',[0.92, 0.1, 0.02, 0.80], 'YDir','normal');
        gradient = linspace(1, length(LF_Shank(end).SortedDepths), length(LF_Shank(end).SortedDepths))';
        imagesc(cbAx, gradient);
        colormap(cbAx, flip(colorMap));  % reuse whatever ‘colorMap’ you built for the deepest shank
        cbAx.XTick = [];
        ylabel(cbAx, 'Depth (µm)', 'FontSize', 12);

        % Save figure
        img_name = sprintf('AllStim_%s_%s_AllShanks_check.png', rat_no, session_name);
        saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
        close;
    end
end

%% Mean IL vs PL per shank, stacked by shank (PL_regions_line)
if PL_regions_line
    for sess = 1:2
        session_name = STIM.epochs{sess};
         validStimIdx = validStimIdx_check{sess};
        nValidStims  = numel(validStimIdx);
        if nValidStims == 0
            warning('Session %s: All stimuli excluded due to noise. Skipping IL vs PL per shank plot.', session_name);
            continue;
        end

        % Build concatenated x-axis offsets
        stim_offsets = (0:(nValidStims-1)) * window_duration;
        % Precompute concatenated x-axis
        x_concat = nan(1, nValidStims * nPts);
        for iVal = 1:nValidStims
            idxs = (iVal-1)*nPts + (1:nPts);
            x_concat(idxs) = x_time + stim_offsets(iVal);
        end

        % Create figure with tiled layout: 4 rows × 1 column
        figure('Position', [50, 50, 2300, 1400],'Visible','off');
        t = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        sgtitle(sprintf('%s Session %s: IL vs PL per Shank', rat_no, session_name), 'FontSize', 20);

        for shank = 1:nShanks
            IL_end   = LF_Shank(shank).IL_end_chan;
            PL_start = IL_end + PRM_IL_PL_gap;

            % Preallocate concatenated mean & std arrays for this shank
            IL_mean_concat = nan(1, nValidStims * nPts);
            IL_std_concat  = nan(1, nValidStims * nPts);
            PL_mean_concat = nan(1, nValidStims * nPts);
            PL_std_concat  = nan(1, nValidStims * nPts);

            for iVal = 1:nValidStims
                iStim = validStimIdx(iVal);
                IL_data = squeeze(LF_Shank(shank).V{sess}(iStim,    1:IL_end,    x_start:x_end));  % [IL_ch x nPts]
                PL_data = squeeze(LF_Shank(shank).V{sess}(iStim, PL_start:end, x_start:x_end));  % [PL_ch x nPts]
                mean_IL = mean(IL_data, 1);
                std_IL  = std(IL_data,  0, 1);
                mean_PL = mean(PL_data, 1);
                std_PL  = std(PL_data,  0, 1);
                idxs = (iVal-1)*nPts + (1:nPts);
                IL_mean_concat(idxs) = mean_IL;
                IL_std_concat(idxs)  = std_IL;
                PL_mean_concat(idxs) = mean_PL;
                PL_std_concat(idxs)  = std_PL;
            end
            % Calc y lim
            ymin2 = min([IL_mean_concat - IL_std_concat, PL_mean_concat - PL_std_concat]);
            ymax2 = max([IL_mean_concat + IL_std_concat, PL_mean_concat + PL_std_concat]);
            ylim([ymin2, ymax2]);

            nexttile;
            hold on;
            ax = gca;
            title(sprintf('Shank %d', shank), 'FontSize', 16);
            if shank == nShanks
                xlabel('Stimulation Amplitude (\muA)', 'FontSize', 14);
                xticks(stim_offsets);
                xticklabels(arrayfun(@(c) sprintf('%.0f', stim_currents(c)), validStimIdx, 'UniformOutput', false));
            else
                ax.XTick = [];
            end
            ylabel('Voltage (mV)', 'FontSize', 14);
           
            set(gca, 'FontSize', 12);

            % Shade IL region (mean ± STD)
            x_patch    = [x_concat, fliplr(x_concat)];
            y_patch_IL = [IL_mean_concat + IL_std_concat, fliplr(IL_mean_concat - IL_std_concat)];
            patch(x_patch, y_patch_IL, IL_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

            % Shade PL region (mean ± STD)
            y_patch_PL = [PL_mean_concat + PL_std_concat, fliplr(PL_mean_concat - PL_std_concat)];
            patch(x_patch, y_patch_PL, PL_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

            % Plot mean IL and PL
            plot(x_concat, IL_mean_concat, 'Color', IL_col, 'LineWidth', 2);
            plot(x_concat, PL_mean_concat, 'Color', PL_col, 'LineWidth', 2);

            % Plot relative 0‐stim line (red) for each stim
            yl = ylim;
            for iVal = 1:nValidStims
                offset = stim_offsets(iVal);
                plot([offset offset], yl, 'Color', 'r', 'LineWidth', 1.5);
            end

            % Plot black vertical lines delimiting each stim window
            for iVal = 1:nValidStims
                start_line = stim_offsets(iVal) + x_time(1);
                end_line   = stim_offsets(iVal) + x_time(end);
                plot([start_line start_line], yl, 'Color', 'k', 'LineWidth', 1);
                plot([end_line   end_line],   yl, 'Color', 'k', 'LineWidth', 1);
            end

            % Legend only on bottom tile
            if shank == nShanks
                legend({'IL ± STD','PL ± STD','Infralimbic','Prelimbic','Zero‐Stim'}, 'FontSize', 12, 'Location', 'northeast');
            end

            hold off;
        end

        % Save figure
        img_name = sprintf('IL_PL_PerShank_%s_%s_AllShanks_check.png', rat_no, session_name);
        saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
        close;
    end
end

%% Mean IL vs PL averaged across shanks (PL_region_avg)
if PL_region_avg
    for sess = 1:2
        session_name = STIM.epochs{sess};
         validStimIdx = validStimIdx_check{sess};
        nValidStims  = numel(validStimIdx);
        if nValidStims == 0
            warning('Session %s: All stimuli excluded due to noise. Skipping regional average plot.', session_name);
            continue;
        end

        stim_offsets = (0:(nValidStims-1)) * window_duration;
        % Preallocate matrices to collect per‐shank concatenated IL/PL
        IL_all = [];
        PL_all = [];

        for shank = 1:nShanks
            IL_end   = LF_Shank(shank).IL_end_chan;
            PL_start = IL_end + PRM_IL_PL_gap;
            IL_concat = nan(1, nValidStims * nPts);
            PL_concat = nan(1, nValidStims * nPts);

            for iVal = 1:nValidStims
                iStim = validStimIdx(iVal);
                IL_data = squeeze(LF_Shank(shank).V{sess}(iStim,    1:IL_end,    x_start:x_end));
                PL_data = squeeze(LF_Shank(shank).V{sess}(iStim, PL_start:end, x_start:x_end));
                IL_concat((iVal-1)*nPts + (1:nPts)) = mean(IL_data, 1);
                PL_concat((iVal-1)*nPts + (1:nPts)) = mean(PL_data, 1);
            end
            IL_all = [IL_all; IL_concat];
            PL_all = [PL_all; PL_concat];
        end

        IL_grand = mean(IL_all, 1);
        PL_grand = mean(PL_all, 1);
        IL_std   = std(IL_all,  0, 1);
        PL_std   = std(PL_all,  0, 1);
        maxValSession = max([IL_grand + IL_std, IL_grand - IL_std, PL_grand + PL_std, PL_grand - PL_std]);
        %ylim
        ymin3 = min([IL_grand - IL_std, PL_grand - PL_std]);
        ymax3 = max([IL_grand + IL_std, PL_grand + PL_std]);
        ylim([ymin3, ymax3]);
        % Build concatenated x-axis
        x_concat = nan(1, nValidStims * nPts);
        for iVal = 1:nValidStims
            idxs = (iVal-1)*nPts + (1:nPts);
            x_concat(idxs) = x_time + stim_offsets(iVal);
        end

        % Create figure for regional average
        figure('Position', [50, 50, 2400, 1800],'Visible','off');
        hold on;
        title(sprintf('%s Session %s: Average IL vs PL Across Shanks', rat_no, session_name), 'FontSize', 20);
        xlabel('Stimulation Amplitude (\muA)', 'FontSize', 16);
        ylabel('Voltage (mV)', 'FontSize', 16);
        set(gca, 'FontSize', 14);
        % Shade IL region
        x_patch    = [x_concat, fliplr(x_concat)];
        y_patch_IL = [IL_grand + IL_std, fliplr(IL_grand - IL_std)];
        patch(x_patch, y_patch_IL, IL_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        % Shade PL region
        y_patch_PL = [PL_grand + PL_std, fliplr(PL_grand - PL_std)];
        patch(x_patch, y_patch_PL, PL_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        % Plot mean IL and PL
        plot(x_concat, IL_grand, 'Color', IL_col, 'LineWidth', 3);
        plot(x_concat, PL_grand, 'Color', PL_col, 'LineWidth', 3);

        % Plot relative 0‐stim line (red) for each stim
        yl = ylim;
        for iVal = 1:nValidStims
            offset = stim_offsets(iVal);
            plot([offset offset], yl, 'Color', 'r', 'LineWidth', 2);
        end

        % Plot black vertical lines delimiting each stim window
        for iVal = 1:nValidStims
            start_line = stim_offsets(iVal) + x_time(1);
            end_line   = stim_offsets(iVal) + x_time(end);
            plot([start_line start_line], yl, 'Color', 'k', 'LineWidth', 1.5);
            plot([end_line   end_line],   yl, 'Color', 'k', 'LineWidth', 1.5);
        end

        xticks(stim_offsets);
        xticklabels(arrayfun(@(c) sprintf('%.0f', stim_currents(c)), validStimIdx, 'UniformOutput', false));
        xlim([x_concat(1), x_concat(end)]);
        legend({'IL \pm STD','PL \pm STD','Infralimbic','Prelimbic','Zero‐Stim'}, 'FontSize', 14, 'Location', 'northeast');

        img_name = sprintf('Avg_IL_PL_AcrossShanks_%s_%s_check.png', rat_no, session_name);
        saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
        close;
    end
end

end
