function[]= NP2_Plot_LF(rat_no,data_filepath,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function[]= NP2_Plot_LF (rat_no,data_filepath,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT: Needs filepath to rat_DATA.mat file

% The function loads plots NP2 data in the following conditions
% - Heat map of  all channels of raw data for every 5th stimulation condition
% - Heat map of  all channels with the average of every 5 consecutive stim conditions
% - Voltage Trace of Average and Stdev of all channels for every 5th stimulation condition
% - Voltage Trace Average and Stdev of  all channels with the average of every 5 consecutive stim conditions
% - IL-PL Separated Voltage Trace of Average and Stdev of all channels for every 5th stimulation condition
% - IL-PL Separated Voltage Trace Average and Stdev of  all channels with the average of every 5 consecutive stim conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS TO CHANGE
% The default runs all plots
PL_heatmap_raw=true;
PL_heatmap_avg = false;
PL_line_trace =true;
PL_line_trace_avg = false;
PL_mPFC_Regions =true;
PL_mPFC_Regions_avg = false;

%Parameters on data set
PRM_IMG_DIR='C:\DATA\Neuropixel\NP2\Processed\Fig_New';
plot_preStim=-50; %100ms
plot_postStim=200; % 400ms
Extract_varargin; % overrides the defaults specified above.

%% Initialize
load (data_filepath, 'PRM' ,'STIM' ,'LF_Shank');
time_start=-1*PRM.time_extract_preStim*1000;  %Convert to ms for the figure
time_end=PRM.time_extract_postStim*1000;
x_time_all = linspace(time_start,time_end, size(LF_Shank(1).V{1,1}, 3));
stim_currents=STIM.sorted_current_uA;

%Get indices of plot=prestim times
[~,x_start]=min(abs(x_time_all-plot_preStim));
[~,x_end]=min(abs(x_time_all-plot_postStim));
x_time=x_time_all(x_start:x_end);

%% Raw Heat Map
if PL_heatmap_raw
    for sess = 1:2
        for stim = 1:length(stim_currents)  % Loop through every 5th stimulus

            % Plotting Raw Traces for Each Shank
            figure;
            set(gcf, 'Position', [100, 100, 1200, 600]);
            t2 = tiledlayout(1, 4, "TileSpacing", "compact", "Padding", "compact");
            caxisValues = [];

            for shank = 1:PRM.shank_no
                ax=nexttile(shank);
                %imagesc(x_time, LF_Shank(shank).SortedDepths, squeeze(LF_Shank(shank).V{sess}(stim, :, :)));
                imagesc(squeeze(LF_Shank(shank).V{sess}(stim, :, x_start:x_end)));
                set(gca, 'YDir', 'normal');
                xlabel('Time (ms)','FontSize',13);
                title(sprintf('Shank %d', shank),'FontSize',14);

                %Setting axis
                xt=round(linspace(x_time(1), x_time(end),numel(ax.XTick)+1));
                yt=flip(round(linspace(LF_Shank(shank).SortedDepths(end),LF_Shank(shank).SortedDepths(1),numel(ax.YTick))/100)*100);
                set(gca,'XTickLabel',xt(2:end),'FontSize',12)
                % Draw Vertical Line at Stim Time/ IL PL Line
                [~, xIndex] = min(abs(x_time));
                line([xIndex xIndex], ylim, 'Color', '#A2142F', 'LineWidth', 2);
                yline(LF_Shank(shank).IL_end_chan, '--k', 'LineWidth', 2);

                if shank > 1
                    set(gca,  'YTickLabel', []);%Disable ytick labels for all but first plot
                else
                    set(gca,'YTickLabel',yt,'FontSize',12);
                    ylabel('Depth (um)','FontSize',14);
                end
                caxisValues = [caxisValues; clim(ax)];  % Collect caxis values from each subplot
            end
            %Color Bar
            for ax = t2.Children'
                clim(ax, [min(caxisValues(:,1)), max(caxisValues(:,2))]);
                %clim(ax, [-5, 5]);
            end

            % Add a single colorbar to the entire layout
            cb = colorbar(ax, 'eastoutside');
            cb.Layout.Tile = 'east';  % This places the colorbar outside the east tile of the layout
            cb.Label.String = 'Response Magnitude (mV)';
            cb.Label.FontSize=14;

            sgtitle(sprintf('%s mPFC Response to %s-Stim:%.0fuA', rat_no, STIM.epochs{sess},stim_currents(stim)),'FontSize',16);
            img_name = sprintf('%s_%s_%.0f_Raw.png', rat_no, STIM.epochs{sess}, stim_currents(stim));
            saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
            close;
        end
    end
end

% %% Average Heat Map
% if PL_heatmap_avg
%     for sess = 1:2
%
%         for stim =1 :PRM.avg_no
%             % Plotting Averages for Each Shank
%             figure;
%             set(gcf, 'Position', [100, 100, 1200, 600]);
%             t2 = tiledlayout(1, 4, "TileSpacing", "compact", "Padding", "compact");
%             caxisValues = [];
%             for shank = 1:PRM.shank_no
%                 ax=nexttile(shank);
%                 imagesc(squeeze(LF_Shank(shank).Avg{sess}(stim, :, :)));
%                 set(gca, 'YDir', 'normal');
%                 xlabel('Time (ms)','FontSize',13);
%                 title(sprintf('Shank %d', shank),'FontSize',14);
%
%                 %Setting axis
%                 xt=linspace(x_time(1), x_time(end),numel(ax.XTick));
%                 yt=flip(round(linspace(LF_Shank(shank).SortedDepths(end),LF_Shank(shank).SortedDepths(1),numel(ax.YTick))/100)*100);
%                 set(gca,'XTickLabel',xt,'FontSize',12)
%                 % Draw Vertical Line at Stim Time & IL PL Line
%                 [~, xIndex] = min(abs(x_time));
%                 line([xIndex xIndex], ylim, 'Color', '#A2142F', 'LineWidth', 2);
%                 yline(LF_Shank(shank).IL_end_chan, '--k', 'LineWidth', 2);
%
%                 if shank > 1
%                     set(gca,  'YTickLabel', []);%Disable ytick labels for all but first plot
%                 else
%                     set(gca,'YTickLabel',yt,'FontSize',12);
%                     ylabel('Depth (um)','FontSize',14);
%                 end
%                 caxisValues = [caxisValues; clim(ax)];  % Collect caxis values from each subplot
%             end
%             %Color Bar
%             for ax = t2.Children'
%                 clim(ax, [min(caxisValues(:,1)), max(caxisValues(:,2))]);
%             end
%
%             % Add a single colorbar to the entire layout
%             cb = colorbar(ax, 'eastoutside');
%             cb.Layout.Tile = 'east';  % This places the colorbar outside the east tile of the layout
%             cb.Label.String = 'Response Magnitude (mV)';
%             cb.Label.FontSize=14;
%
%             sgtitle(sprintf('%s mPFC Average Response to %s-Stim:%.0f-%.0fuA', rat_no,STIM.epochs{sess},...
%                 stim_currents(stim*5-4),stim_currents(stim*5)),'FontSize',16);
%             img_name = sprintf('%s_%s_%.0f_Average.png', rat_no, STIM.epochs{sess}, stim_currents(stim));
%             saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
%             close;
%
%         end
%     end
% end

%% Line Trace

if PL_line_trace
    % Generate a color map
    colors = [
        0.35, 0.05, 0.7     % Purple
        0.25, 0.41, 0.88; % Blue
        0.15, 0.5, 0.6 %
        0.55    0.7    0.4 %Light green
        0.9    0.85    0.26 %yellow
        0.88, 0.35, 0.11]; %Orange

    % Determine the number of channels across all shanks
    max_nChan = max([length(LF_Shank(1).SortedDepths),length(LF_Shank(2).SortedDepths),length(LF_Shank(3).SortedDepths),length(LF_Shank(4).SortedDepths)]);
    % Generate a piecewise linear colormap
    num_colors = size(colors, 1);
    x = linspace(0, 1, num_colors);  % Original positions
    xq = linspace(0, 1, max_nChan);  % Interpolated positions

    % Interpolate each RGB channel
    colorMap = zeros(max_nChan, 3);
    for i = 1:3
        colorMap(:, i) = interp1(x, colors(:, i), xq, 'linear');
    end

    for sess = 1:2
        for stim = 1:length(stim_currents)  % Loop through every 5th stimulus
            figure;
            set(gcf,'Position',[100,100,1250,600])
            t2=tiledlayout(1, 4, "TileSpacing", "compact", "Padding", "compact");
            set(t2, 'InnerPosition', [0.06 0.1 0.85 .80]);
            axHandles = gobjects(1, PRM.shank_no);
            for shank = 1:PRM.shank_no
                axHandles(shank) =nexttile(shank);
                hold on;
                % Plot each channel with color from the color map
                for chan = 1:length(LF_Shank(shank).SortedDepths)
                    dt = squeeze(LF_Shank(shank).V{sess}(stim, chan, x_start:x_end));
                    plot(x_time, dt, 'LineWidth', 1, 'Color', colorMap(chan, :));
                end
                set(gca,'FontSize',12);
                xlim([min(x_time), max(x_time)]);
                xlabel('Time (ms)','FontSize',14);
                title(sprintf('Shank %d', shank));

            end
            linkaxes(axHandles, 'y');
            ylim([-3 3]);
            for shank = 1:PRM.shank_no
                axHandles(shank) = nexttile(shank);
                hold on;

                % --- compute vertical spacing so traces don't overlap ---
                nChan = length(LF_Shank(shank).SortedDepths);
                % pull out all data for this session, this shank, across time & stims
                dataWindow = reshape(LF_Shank(shank).V{sess}(:,:, x_start:x_end), [], 1);
                dataRange  = max(dataWindow) - min(dataWindow);
                offset_spacing = 0.075;%dataRange * 1.2;  % 20% extra gap

                % --- plot each channel with its own offset ---
                for chan = 1:nChan
                    dt = squeeze(LF_Shank(shank).V{sess}(stim, chan, x_start:x_end));
                    dt_offset = dt + (chan-1) * offset_spacing;
                    plot(x_time, dt_offset, 'LineWidth', 1, 'Color', colorMap(chan, :));
                end

                % --- label y‑axis in units of depth (or channel index) ---
                yTicks   = (0:(nChan-1)) * offset_spacing;
                yLabels  = LF_Shank(shank).SortedDepths;    % in μm; or use 1:nChan if you prefer
                set(gca, 'YTick', yTicks, 'YTickLabel', []);

                xlim([min(x_time), max(x_time)]);
                xlabel('Time (ms)','FontSize',14);
                title(sprintf('Shank %d', shank),'FontSize',14);
            end

            % link y‑limits so all shanks share the same vertical span (optional)
            linkaxes(axHandles, 'y');

            for shank=1:PRM.shank_no
                nexttile(shank)
                line([min(abs(x_time_all)) min(abs(x_time_all))], ylim, 'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time
                if shank ==1
                    ylabel('Voltage(mV)', 'FontSize', 14);
                end
            end
            %Colorbar
            cbAx = axes('Position', [0.915 0.1 0.02 0.80], 'YDir', 'normal'); % Position [left bottom width height]
            gradient = linspace(1, length(LF_Shank(4).SortedDepths), length(LF_Shank(4).SortedDepths))';
            imagesc(cbAx, gradient);
            colormap(cbAx, flip(colorMap));
            yt=(round(linspace(LF_Shank(shank).SortedDepths(end),LF_Shank(shank).SortedDepths(1),numel(cbAx.YTick))/100)*100);
            set(gca,'YTickLabel',yt,'FontSize',10,'YAxisLocation','right','YTickLabelRotation',30);
            ylabel('Depth (um)','FontSize',14);
            % Remove x-ticks and x-label from the color bar
            xticks([]);
            xlabel(cbAx, '');
            %Main Title
            sgtitle(sprintf('%s mPFC Voltage Trace to %s-Stim:%.0fuA', rat_no, STIM.epochs{sess},stim_currents(stim)),...
                'FontSize',16);
            pubify_figure_axis_robust(16,16);
            img_name = sprintf('ANEW%s_%s_%.0f_Line.png', rat_no, STIM.epochs{sess}, stim_currents(stim));
            saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
            close;
        end
    end
end
%% Average Line Trace

if PL_line_trace_avg
    lightBlue = [0.500, 0.700, 0.770];  % Light blue color
    darkBlue = [0.1, 0.1, 0.3];  % Dark blue color
    for sess = 1:2
        for stim = 1:PRM.avg_no
            figure;
            set(gcf,'Position',[100,100,1250,600])
            t2=tiledlayout(1, 4, "TileSpacing", "compact", "Padding", "compact");
            axHandles = gobjects(1, PRM.shank_no);
            set(t2, 'InnerPosition', [0.06 0.1 0.85 .80]);
            for shank = 1:PRM.shank_no
                axHandles(shank) =nexttile(shank);
                hold on;
                colorMap = [linspace(darkBlue(1), lightBlue(1), length(LF_Shank(shank).SortedDepths))', ...
                    linspace(darkBlue(2), lightBlue(2), length(LF_Shank(shank).SortedDepths))', ...
                    linspace(darkBlue(3), lightBlue(3), length(LF_Shank(shank).SortedDepths))'];

                % Plot each channel with color from the color map
                for chan = 1:length(LF_Shank(shank).SortedDepths)
                    dt = squeeze(LF_Shank(shank).Avg{sess}(stim, chan, :));
                    plot(x_time, dt, 'LineWidth', 0.5, 'Color', colorMap(chan, :));
                end
                set(gca,'FontSize',12);
                xlim([min(x_time), max(x_time)]);
                line([min(abs(x_time)) min(abs(x_time))], ylim, 'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time
                xlabel('Time (ms)','FontSize',14);
                title(sprintf('Shank %d', shank));

            end
            %Y axis and plots
            linkaxes(axHandles, 'y');
            for shank=1:PRM.shank_no
                nexttile(shank)
                line([min(abs(x_time)) min(abs(x_time))], ylim, 'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time
                if shank ==1
                    ylabel('Voltage(mV)', 'FontSize', 14);
                else
                    set(gca,  'YTickLabel', []);
                end
            end
            %Color Bar
            cbAx = axes('Position', [0.915 0.1 0.02 0.80], 'YDir', 'normal'); % Position [left bottom width height]
            gradient = linspace(1, length(LF_Shank(4).SortedDepths), length(LF_Shank(4).SortedDepths))';
            imagesc(cbAx, gradient);
            colormap(cbAx, [linspace(lightBlue(1), darkBlue(1), length(gradient))', ...
                linspace(lightBlue(2), darkBlue(2), length(gradient))', ...
                linspace(lightBlue(3), darkBlue(3), length(gradient))']);

            yt=(round(linspace(LF_Shank(shank).SortedDepths(end),LF_Shank(shank).SortedDepths(1),numel(cbAx.YTick))/100)*100);
            set(gca,'YTickLabel',yt,'FontSize',10,'YAxisLocation','right','YTickLabelRotation',30);
            ylabel('Depth (um)','FontSize',14);
            % Remove x-ticks and x-label from the color bar
            xticks([]);
            xlabel(cbAx, '');
            %Main Title
            sgtitle(sprintf('%s mPFC Average Voltage Trace to %s-Stim:%.0f-%.0fuA', rat_no,STIM.epochs{sess},...
                stim_currents(stim*5-4),stim_currents(stim*5)),'FontSize',16);
            img_name = sprintf('%s_%s_%.0f_Line_Avg.png', rat_no, STIM.epochs{sess}, stim_currents(stim));
            saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
            close;
        end
    end
end
%% mPFC -Infralimbic Prelimbic Separation
if PL_mPFC_Regions
    PRM_IL_PL_gap=4;
    PL_col = [0.88, 0.35, 0.11]; % Orange
    IL_col = [0.15, 0.5, 0.6]; % Teal
    for sess = 1:2
        for stim = 1:5:length(stim_currents)  % Loop through every 5th stimulus
            % Plotting Raw Traces for Each Shank
            figure;
            set(gcf,'Position',[100,100,1250,600])
            t2=tiledlayout(1, 4, "TileSpacing", "compact", "Padding", "compact");
            axHandles = gobjects(1, PRM.shank_no);
            set(t2, 'InnerPosition', [0.06 0.1 0.85 .80]);
            for shank = 1:PRM.shank_no
                axHandles(shank) =nexttile(shank);
                mean_IL = mean(squeeze(LF_Shank(shank).V{sess}(stim, 1:LF_Shank(shank).IL_end_chan,x_start :x_end)), 1);
                mean_PL = mean(squeeze(LF_Shank(shank).V{sess}(stim, LF_Shank(shank).IL_end_chan+PRM_IL_PL_gap:end, x_start :x_end)), 1);
                hold on;
                h1=plot(x_time,mean_IL,'LineWidth',4,'Color',IL_col);
                h2=plot(x_time,mean_PL,'LineWidth',4,'Color',PL_col);
                %Plot channels

                % Channel plots - vectorized approach
                IL_channels = squeeze(LF_Shank(shank).V{sess}(stim, 1:LF_Shank(shank).IL_end_chan, x_start :x_end));
                PL_channels = squeeze(LF_Shank(shank).V{sess}(stim, LF_Shank(shank).IL_end_chan+PRM_IL_PL_gap:end, x_start :x_end));
                plot(x_time, IL_channels, 'LineWidth', 1.5, 'Color', [IL_col, 0.3]);
                plot(x_time, PL_channels, 'LineWidth', 1.5, 'Color', [PL_col, 0.3]);
                %Axis
                set(gca,'FontSize',12);
                xlim([min(x_time), max(x_time)]);
                line([0, 0], ylim, 'Color', '#A2142F', 'LineWidth', 4);  % Vertical line at x_time = 0
                xlabel('Time (ms)','FontSize',14);
                title(sprintf('Shank %d', shank),'FontSize',14);

                hold off;
            end
            %Y axis link
            linkaxes(axHandles, 'y');
            for shank=1:PRM.shank_no
                nexttile(shank)
                line([min(abs(x_time)) min(abs(x_time))], ylim, 'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time
                if shank ==1
                    ylabel('Voltage(mV)', 'FontSize', 14);
                else
                    set(gca,  'YTickLabel', []);
                end
            end
            %Plot legend properties
            legend([h1, h2], {'Infralimbic', 'Prelimbic'}, 'Location', 'northeast', 'Orientation', 'vertical');
            sgtitle(sprintf('%s Infralimbic vs Prelimbic Response to %s-Stim:%.0f', rat_no, STIM.epochs{sess},stim_currents(stim)),'FontSize',16);
            img_name = sprintf('IL_PL_%s_%s_%.0f_Line.png', rat_no, STIM.epochs{sess}, stim_currents(stim));
            saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
            close;
        end
    end
end

%% Average IL-PL
% if PL_mPFC_Regions_avg
%     PRM_IL_PL_gap=4;
%     PL_col = [0.2, 0.5, 0.6];  % Teal
%     IL_col = [0.9, 0.6, 0.2]; % Orange
%     for sess = 1:2
%         for stim = 1:PRM.avg_no
%             % Plotting Raw Traces for Each Shank
%             figure;
%             set(gcf,'Position',[100,100,1250,600])
%             t2=tiledlayout(1, 4, "TileSpacing", "compact", "Padding", "compact");
%             axHandles = gobjects(1, PRM.shank_no);
%             set(t2, 'InnerPosition', [0.06 0.1 0.85 .80]);
%             for shank = 1:PRM.shank_no
%                 axHandles(shank) =nexttile(shank);
%                 hold on;
%                 %Mean
%                 mean_IL = mean(squeeze(LF_Shank(shank).Avg{sess}(stim, 1:LF_Shank(shank).IL_end_chan, :)), 1);
%                 mean_PL = mean(squeeze(LF_Shank(shank).Avg{sess}(stim, LF_Shank(shank).IL_end_chan+PRM_IL_PL_gap:end, :)), 1);
%                 h1=plot(x_time,mean_IL,'LineWidth',4,'Color',IL_col);
%                 h2=plot(x_time,mean_PL,'LineWidth',4,'Color',PL_col);
%                 %STD
%                 std_IL = std(squeeze(LF_Shank(shank).Avg{sess}(stim, 1:LF_Shank(shank).IL_end_chan, :)),0, 1);
%                 std_PL = std(squeeze(LF_Shank(shank).Avg{sess}(stim, LF_Shank(shank).IL_end_chan+PRM_IL_PL_gap:end, :)),0, 1);
%                 fill([x_time, fliplr(x_time)],[mean_IL + std_IL, fliplr(mean_IL - std_IL)],IL_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%                 fill([x_time, fliplr(x_time)], [mean_PL + std_PL,fliplr(mean_PL - std_PL)], PL_col,'FaceAlpha', 0.3, 'EdgeColor', 'none');
%                 %Axis
%                 set(gca,'FontSize',12);
%                 xlim([min(x_time), max(x_time)]);
%                 line([0, 0], ylim, 'Color', '#A2142F', 'LineWidth', 4);  % Vertical line at x_time = 0
%                 xlabel('Time (ms)','FontSize',14);
%                 title(sprintf('Shank %d', shank),'FontSize',14);
%                 hold off;
%             end
%              %Y axis link
%             linkaxes(axHandles, 'y');
%             for shank=1:PRM.shank_no
%                 nexttile(shank)
%                 line([min(abs(x_time)) min(abs(x_time))], ylim, 'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time
%                 if shank ==1
%                     ylabel('Voltage(mV)', 'FontSize', 14);
%                 else
%                     set(gca,  'YTickLabel', []);
%                 end
%             end
%             legend([h1, h2], {'Infralimbic', 'Prelimbic'}, 'Location', 'northeast', 'Orientation', 'vertical');
%             sgtitle(sprintf('%s Infralimbic vs Prelimbic Response to %s-Stim:%.0f-%.0f', rat_no,STIM.epochs{sess}, ...
%                 stim_currents(stim*5-4),stim_currents(stim*5)),'FontSize',16);
%             img_name = sprintf('IL_PL_%s_%s_%.0f_Line_Avg.png', rat_no, STIM.epochs{sess}, stim_currents(stim));
%             saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
%             close;
%         end
%     end
% end
%%
if PL_mPFC_Regions_avg
    PRM_IL_PL_gap = 4;
    IL_col = [0.35, 0.05, 0.7];
    PL_col = [0.88, 0.35, 0.11];
    for sess = 1:2
        % Plotting Raw Traces for all Shanks on the same graph
        figure;
        set(gcf, 'Position', [100, 100, 1250, 600])
        hold on;

        for shank = 1:4
            % Calculate average across all stim conditions for IL and PL
            avg_IL = mean(mean(LF_Shank(shank).Avg{sess}(:, 1:LF_Shank(shank).IL_end_chan, :), 1), 2); % Mean across dims 1 and 2
            avg_PL = mean(mean(LF_Shank(shank).Avg{sess}(:, LF_Shank(shank).IL_end_chan + PRM_IL_PL_gap:end, :), 1), 2); % Mean across dims 1 and 2

            % Squeeze to remove singleton dimensions for plotting
            avg_IL = squeeze(avg_IL);
            avg_PL = squeeze(avg_PL);

            % Mean plotting
            h1 = plot(x_time, avg_IL, 'LineWidth', 4, 'Color', IL_col);
            h2 = plot(x_time, avg_PL, 'LineWidth', 4, 'Color', PL_col);

            % Standard Deviation across dimensions 1 and 2
            std_IL = squeeze(std(mean(LF_Shank(shank).Avg{sess}(:, 1:LF_Shank(shank).IL_end_chan, :), 1), 0, 2));
            std_PL = squeeze(std(mean(LF_Shank(shank).Avg{sess}(:, LF_Shank(shank).IL_end_chan + PRM_IL_PL_gap:end, :), 1), 0, 2));

            % Fill plots for standard deviations
            fill([x_time, fliplr(x_time)], [avg_IL + std_IL; flipud(avg_IL - std_IL)], IL_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            fill([x_time, fliplr(x_time)], [avg_PL + std_PL; flipud(avg_PL - std_PL)], PL_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        end
        % Axis settings
        set(gca, 'FontSize', 12);
        xlim([-50, 150]);
        ylim([-2.5 1.5 ])
        line([0, 0], ylim, 'Color', '#A2142F', 'LineWidth', 4);  % Vertical line at x_time = 0
        xlabel('Time (ms)', 'FontSize', 14);
        ylabel('Voltage (mV)', 'FontSize', 14);
        %title(sprintf('Shank %d', shank), 'FontSize', 14);
        hold off;

        % Legend and title
        legend([h1, h2], {'Infralimbic', 'Prelimbic'}, 'Location', 'northeast', 'Orientation', 'vertical');
        sgtitle(sprintf('Old Rats Average Response to vHC Stim'), 'FontSize', 16);

        % Save image
        img_name = sprintf('IL_PL_%s_%s_Avg_AllStim_Shank%d.png', rat_no, STIM.epochs{sess}, shank);
        saveas(gcf, fullfile(PRM_IMG_DIR, img_name));
        close;
    end
end

