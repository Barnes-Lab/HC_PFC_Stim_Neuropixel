function[]= NP2_LF_DATA_LOAD_FROM_AP (rat,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function NP2_LF_DATA_LOAD(PRM.RAT_NO,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function loads Neuropixel Data for vHC_STIM Project for each ratbn
% -Loads Stim params from stim file name
% -Loads Chan map -name is kilosortChanMap
% -Loads Meta Data
% -Loads LF data post stim for all channels across each condition - assumes
%   iHC, vHC, Burst Stim Conditions (Time is 150 ms can change parameter)
% -Removes the Noise by subtracting data from the time before stim
% -Sorts LF data on the basis of shank and gets the depth for each site
% -Saves data in .mat file for further processing
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS TO CHANGE
%SOFTWARE DIRS
PRM_TPRIME_DIR='C:\DATA\Neuropixel\Software\TPrime-win';

%DATA DIRS
PRM_DATA_DIR='D:\DATA\Neuropixel\NP2';
PRM_STORAGE_DIR='D:\DATA\Neuropixel\NP2\Processed'; %Where data is stored and worked form later
PRM_SUB_FOLDER = {'mPFC_L23_L5_NP2_g0'};
PRM_SUB_FOLDER_catgt={'catgt_mPFC_L23_L5_NP2_g0'};
PRM_STIM_CHANNEL=5;% Default channel on which the stim currents are extracted for TPrime/etc

% Data Parameters
PRM.n_stim=25; %no. of stims per burst
PRM.epochs={'iHC','vHC','BurstStim'}; %Different epochs of experiment
PRM.chanFile='mPFC_L23_L5_NP2_g0_tcat.imec0.ap_kilosortChanMap.mat';
PRM.time_extract_postStim=.6; %Time in sec to extract post stim 400ms
PRM.time_extract_preStim=1; %Time in sec to extract post stim 200ms
PRM.pre_plot=0.2;
PRM.post_plot=0.4 ;% 400 ms

%SHANK AND AVERAGING DATA SPECS
PRM.shank_no=4; %No. of shanks
PRM.avg_no=5; %Averages across 5 sessions
PRM.max_depth=5125; %Depth of probe In micrometers (Calculated from sterotax coordinates -tip)
%IL PL Separation
PRM.PL_depth=3600; %Depth of PL-IL transition in micrometers

Extract_varargin; % overrides the defaults specified above.
%% STIM PARAMS LOAD
%Fig dirs
PRM.RAT_NO=rat;
PRM.RAT_DATA_DIR=fullfile(PRM_DATA_DIR,PRM.RAT_NO);
PRM.RAT_STORAGE_DIR=fullfile(PRM_STORAGE_DIR,PRM.RAT_NO);
PRM.IMG_DIR=fullfile(PRM.RAT_STORAGE_DIR,'Figures');

if ~exist(PRM.IMG_DIR,"dir")
    mkdir(PRM.IMG_DIR)
end

%Load params from Original Stim File
params_path=dir(fullfile(PRM.RAT_DATA_DIR,'params*.mat'));
STIM_PRM=load(fullfile(params_path.folder,params_path.name));
stim_currents=STIM_PRM.PRM.output_current_array;
[STIM.sorted_current,STIM.sorted_idx] = sort(STIM_PRM.PRM.output_current_array);
STIM.sorted_current_uA=(STIM.sorted_current-0.0986)./0.0113;                            %%% CHANGE THIS VALUE IF USING A DIFF STIMULATION



%% GET SESS INFO, RUN TPRIME
obj1 = SGLX_Class;
%GET SESSION FILES AND RUN TPRIME
SESS.ROOT=fullfile(PRM.RAT_DATA_DIR,PRM_SUB_FOLDER); 
SESS.ROOT_catgT=fullfile(PRM.RAT_DATA_DIR,PRM_SUB_FOLDER_catgt);                   %Root session folder
SESS.FILE_DIR=fullfile(SESS.ROOT,strrep(PRM_SUB_FOLDER, 'g0','g0_imec0'));          %Sub folder with .bin files
SESS.LF_meta_fname=strrep(PRM_SUB_FOLDER,'g0','g0_tcat.imec0.ap.meta');%LF meta name
SESS.LF_meta=obj1.ReadMeta(SESS.LF_meta_fname{1},SESS.ROOT_catgT{1});             %meta data


%Run TPrime
% NPXL_Sync_With_TPrime_SS('PRM_TPRIME_DIR',PRM_TPRIME_DIR,'PRM_ROOT_DATA_DIR', SESS.ROOT{1},... 
%     'PRM_evt_files_extension',{sprintf('*nidq.x*_%d_0.txt',PRM_STIM_CHANNEL)},'PRM_EVT_generate',false);

%Get Stim Times
STIM.fname=fullfile(SESS.ROOT,sprintf('synced_%s_tcat.nidq.xa_%d_0.txt',PRM_SUB_FOLDER{1},PRM_STIM_CHANNEL));
stimes=readmatrix(STIM.fname{1});
%SPLIT INTO EPOCHS
STIM.epochs=PRM.epochs; %Change if order was switched
STIM.times{1}=stimes(1:PRM.n_stim); %iHC
STIM.times{2}=stimes(PRM.n_stim+1:2*PRM.n_stim); %vHC
STIM.times{3}=stimes(2*PRM.n_stim+5:5:end); %Burst Stim vHC

%Load Channel Map Info
ChanMap=load(fullfile(SESS.FILE_DIR{1},PRM.chanFile));
%% Load LFP DATA
%Initiate arrays
ALL_DATA = {};

% LFP session metadata
LFP.meta        = SESS.LF_meta;
LFP.meta_fname  = SESS.LF_meta_fname{1};
LFP.FILE_DIR    = SESS.ROOT_catgT{1};

% Sampling parameters
time_diff = PRM.time_extract_preStim + PRM.time_extract_postStim; % total window per stim (s)
Fs        = str2double(LFP.meta.imSampRate);           % original sampling rate
Fs_ds     = floor(Fs/15);                                % downsampled rate
LFP.nChan = str2double(LFP.meta.nSavedChans);
LFP.nSamp = floor(time_diff * Fs);                      % samples per stim window

ALL_DATA.badChans=[237,299,49];


% Convert to Volts
gain  = str2double(LFP.meta.imChan0apGain);
imax  = str2double(LFP.meta.imMaxInt);
vmax  = str2double(LFP.meta.imAiRangeMax);


% Loop over sessions
for sess = 1:2
    numStim = length(STIM.sorted_current_uA);
    ALL_DATA(sess).Time = nan(numStim, floor(LFP.nSamp/15));
    ALL_DATA(sess).FiltData = nan(numStim, LFP.nChan, floor(LFP.nSamp/15));

    stim_times = STIM.times{sess};
    LFP.binName = strrep(LFP.meta_fname, '.ap.meta', '.ap.bin');

    for stim = 1:length(stim_times)
        % Calculate start sample for this stim segment
        t0    = stim_times(stim) - PRM.time_extract_preStim;
        samp0 = floor(t0 * Fs);

        % Read raw segment [t0 : t0 + time_diff]
        rawSeg = obj1.ReadBin(samp0, LFP.nSamp, LFP.meta, LFP.binName, LFP.FILE_DIR);
        rawSeg_V = rawSeg * vmax / (imax * gain) * 1000;

        %Read DC pre stim 3 secs
        t1 = t0-2; % t0 is already 1 sec prior
        idx0 = floor(t1 * Fs);        % zero-based sample index 
        nWin = floor((t0+1 - t1) * Fs); % number of samples = Fs * 1 s
        % Get data
        chunkData   = obj1.ReadBin(idx0, nWin, LFP.meta, LFP.binName, LFP.FILE_DIR);  % [nChan x nWin]
        chunkV = chunkData * vmax / (imax * gain) * 1000;  % [nChan x nWin]
        dcOffsetMethod1 = mean(chunkV, 2);  %
        clear chunkData chunkV;

        % DC Offset
        data_noDC1 = rawSeg_V - dcOffsetMethod1;  % subtract per-channel offset
        dataFilt1  = downsampleLFP(data_noDC1, Fs, Fs_ds);
        dataSmooth1= SmoothSignal(dataFilt1, Fs_ds);

       
        % Time vector after downsampling
        nDS = size(dataFilt1, 2);
        tvec = t0 + (0:(nDS-1)) / Fs_ds;  % seconds from session start

        ALL_DATA(sess).Time(stim, :) = tvec;
        ALL_DATA(sess).FiltData(stim, :, :) = dataSmooth1;
        
        [~,x_zero]=min(abs(tvec-stim_times(stim)));
        % Clear temporary variables to free memory
        clear rawSeg rawSeg_V data_noDC1 dataFilt1 dataSmooth1 data_noDC2 dataFilt2 dataSmooth2;
        % Plotting if needed to verify 
        % %All channels    
         % if mod(stim,5)==0
         %   figure
         %   hold on;
         %   tiledlayout(1, 4, "TileSpacing", "compact")
         %   %p1
         %   nexttile
         %    hold on
         %   plot(tvec(1900:2400),rawSeg_V(:,28500:15:36000),'k','LineWidth',0.1)
         %   line([tvec(x_zero) tvec(x_zero)], ylim,'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time  %p2   
         %   title('Raw Data')
         %   %p2 
         %   nexttile
         %     hold on
         %   plot(tvec(1900:2400),data_noDC1(:,28500:15:36000),'k','LineWidth',0.1)
         %   line([tvec(x_zero) tvec(x_zero)], ylim,'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time
         %   title('DC subtract')
         %    %p3
         %   nexttile
         %   hold on
         %   plot(tvec(1900:2400),dataFilt1(:,1900:2400),'k','LineWidth',0.1)
         %   line([tvec(x_zero) tvec(x_zero)], ylim,'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time  title('Downsampled')
         %   %p4
         %   nexttile
         %    hold on
         %   plot(tvec(1900:2400),dataSmooth1(:,1900:2400),'k','LineWidth',0.1)
         %  line([tvec(x_zero) tvec(x_zero)], ylim,'Color', '#A2142F', 'LineWidth', 3);% Draw Vertical Line at Stim Time    title('Smoothed')
         %   sgtitle(sprintf('%sFiltered Data %s Stim:%.2f ',PRM.RAT_NO, STIM.epochs{sess}, round(STIM.sorted_current_uA(STIM.sorted_idx(stim)))));
         %   img_name = sprintf('FD_%s_%s_%.0f_.png', PRM.RAT_NO, STIM.epochs{sess}, round(STIM.sorted_current_uA(STIM.sorted_idx(stim))));
         %   saveas(gcf, fullfile(PRM.IMG_DIR, img_name));
         %   close;
         % end
         % % 
    end
end


%% SORTING DATA BY SHANK

% Initialize LF_Shank structure
LF_Shank = struct();
badchans=ALL_DATA.badChans;

%% Processing for Each Shank
for ii = 1:PRM.shank_no
    % Initialize or reset the V and Avg fields for the current shank
    LF_Shank(ii).V = [];
    LF_Shank(ii).Avg = [];
    LF_Shank(ii).SortedDepths=[];

    % Extract shank channels and depths, and remove bad channels
    shank_channels = round(ChanMap.chanMap0ind(ChanMap.kcoords == ii));
    depths = ChanMap.ycoords(ChanMap.kcoords == ii);

    % Identify and remove bad channels for current shank
    valid_idx = ~ismember(shank_channels, badchans);
    valid_shank_channels = shank_channels(valid_idx);
    valid_depths = depths(valid_idx);

    % Sort depths 
    [sorted_depths, depth_sort_idx] = sort(valid_depths); 
    sorted_shank_channels = valid_shank_channels(depth_sort_idx) + 1; % Apply the sorted depths to shank channels (maintains order)
    %Calculates the depth of each tetrode (based on stereotax depth-site depth)
    LF_Shank(ii).SortedDepths =PRM.max_depth- sorted_depths;
    

    %% Processing for Each Session
    for sess = 1:2
        LF_data = ALL_DATA(sess).FiltData;
        % Select LF data based on sorted shank channels
        sel_LF = LF_data(STIM.sorted_idx, sorted_shank_channels, :);
     
        % Store results in LF_Shank
        LF_Shank(ii).V{sess} = sel_LF;
       
        % % Plot Unsorted LFP Data
        % % Define the sample index for the vertical line (e.g., 2000)
        % stim_sample_idx = 2000;
        % time_vector = (0:size(LF_data, 3) - 1) / Fs_ds - PRM.time_extract_preStim;
        % figure('Name', 'Unsorted LFP Data - First Trial', 'Position', [100, 100, 1200, 600]);
        % hold on;
        % for ch = 1:size(LF_data, 2)
        %     plot(time_vector, squeeze(LF_data(3, ch, :)) + ch * 0.1, 'k'); % Offset each channel for clarity
        % end
        % xline(time_vector(stim_sample_idx), 'r', 'LineWidth', 2); % Vertical line at stimulus time
        % title('Unsorted LFP Data - First Trial');
        % xlabel('Time (s)');
        % ylabel('Channel (offset)');
        % hold off;
        % 
        % % Plot Sorted LFP Data
        % figure('Name', 'Sorted LFP Data - First Trial', 'Position', [100, 100, 1200, 600]);
        % hold on;
        % for ch = 1:size(sel_LF, 2)
        %     plot(time_vector, squeeze(sel_LF(1, ch, :)) + ch * 0.1, 'b'); % Offset each channel for clarity
        % end
        % xline(time_vector(stim_sample_idx), 'r', 'LineWidth', 2); % Vertical line at stimulus time
        % title('Sorted LFP Data - First Trial');
        % xlabel('Time (s)');
        % ylabel('Channel (offset)');
        % hold off;

        clear LF_data sel_LF;

    end
end
%Calculate the IL PL Depth difference chaneel
for shank=1:PRM.shank_no
    differences = abs(LF_Shank(shank).SortedDepths - PRM.PL_depth);
    [~, LF_Shank(shank).IL_end_chan] = min(differences);%index of the minimum difference
end    

%% SAVE IN DATA FILE
save_fname=fullfile(PRM.RAT_STORAGE_DIR,sprintf('%s_DATA.mat',PRM.RAT_NO));
save(save_fname,'STIM','ALL_DATA','ChanMap','LF_Shank','PRM');

fprintf('\n Rat Data for %s saved in %s',PRM.RAT_NO,save_fname)

clear ('ALL_DATA','LFP','STIM','LF_Shank','PRM');

end


























