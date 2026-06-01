%% 1) INPUT PARAMETERS 
% FInal Rats

T_PRM.RATS={'10940','10947','10986','10987','10994',...
    '11039','11042','11043','11040','11044'};
T_PRM.AGES={'O','O','Y','Y','O','O','Y','Y','O','Y'};
T_PRM.STIM_CHANNEL={5,5,5,5,4,4,4,4,4}; %Channel on which the stim pulse was sent - need for Tprime

% Data Dirs
T_PRM.TPRIME_DIR='I:\NEUROPIXEL\Software\TPrime-win';
T_PRM.DATA_DIR='E:\NP2'; % Should have the tcat.ap.bin files for each rat
T_PRM.STORAGE_DIR='C:\DATA\NP2\Processed'; %Where data is stored and  worked form later
T_PRM.SPIKE_DIR='C:\DATA\NP2';
T_PRM.SUB_FOLDER = {'mPFC_L23_L5_NP2_g0'};
T_PRM.shank_no=4;

% Make Image Dir
PRM_ALL_IMG_DIR=fullfile(T_PRM.STORAGE_DIR,'LFP_All');
if ~exist(PRM_ALL_IMG_DIR,"dir"); mkdir(PRM_ALL_IMG_DIR); end
T_PRM.IMG_DIR=PRM_ALL_IMG_DIR;

%% 2) Load and Process Data
% This loads the LF data file for the stim and saves it in .mat files
% separated by shank

for rat=1:length(T_PRM.RATS)
    NP2_LF_DATA_LOAD_FROM_AP(T_PRM.RATS{rat},'PRM_TPRIME_DIR',T_PRM.TPRIME_DIR,...
        'PRM_DATA_DIR',T_PRM.DATA_DIR,'PRM_STORAGE_DIR',T_PRM.STORAGE_DIR,...
        'PRM_STIM_CHANNEL',T_PRM.STIM_CHANNEL{rat});
end
%% 2b) %Optional: You can load all the parameter files to check the stim order 
%Additional step is to load and verify with .nidq files (check stim folder)

%Load Stim from data file
for rat = 1:length(T_PRM.RATS)
    data_filepath = fullfile(T_PRM.STORAGE_DIR, T_PRM.RATS{rat}, sprintf('%s_DATA.mat', T_PRM.RATS{rat}));
    load(data_filepath, 'STIM');

    params_path = dir(fullfile(T_PRM.DATA_DIR, T_PRM.RATS{rat}, 'params*.mat'));
    STIM_PRM = load(fullfile(params_path.folder, params_path.name));
    stim_currents = STIM_PRM.PRM.output_current_array;

    % Verify: applying sorted_idx to original must recover sorted_current
    if ~all(stim_currents(STIM.sorted_idx) == STIM.sorted_current)
        fprintf('Check stim order for rat: %s\n', T_PRM.RATS{rat});
    else
        fprintf('Stim order verified for rat: %s\n', T_PRM.RATS{rat});
    end
end

%% 3) Load All Stim Times after TPrime into one file
% This must be run after the Load (2) where the TPrime files are
% saved 
% this output is used for all future analysis

% Combine all stims across rats
ALL_STIM_TIMES = cell(1, length(T_PRM.RATS));
allStimData = [];  % Initialize array for CSV export

for rat = 1:length(T_PRM.RATS)
    data_filepath = fullfile(T_PRM.STORAGE_DIR, T_PRM.RATS{rat}, sprintf('%s_DATA.mat', T_PRM.RATS{rat}));
    load(data_filepath, 'STIM');
    
    % Each session has aligned current and time values
    for s = 1:length(STIM.epochs)-1
        times = STIM.times{s}(:);  
        currents = STIM.sorted_current_uA'; 
        
        [A,idx]=sort(currents);
        currents=currents(idx);
        times=times(idx);
        
        rat_list = repmat(string(T_PRM.RATS{rat}), numel(times), 1);
        session = repmat(string(STIM.epochs{s}), numel(times), 1);
        
       
        tbl = table(rat_list, session, currents(:), times, ...
            'VariableNames', {'Rat', 'Session', 'Current_uA', 'Time_s'});
        
        ALL_STIM_TIMES{rat, s} = tbl;
        allStimData = [allStimData; tbl];
    end
end

% Save .mat and .csv files (USED ALL PLACES DOWNSTREAM)
save(fullfile(T_PRM.STORAGE_DIR, 'ALL_STIM_TIMES_2.mat'), 'ALL_STIM_TIMES');
writetable(allStimData, fullfile(T_PRM.STORAGE_DIR, 'ALL_STIM_TIMES_2.csv'));

fprintf('\nAll Stim times completed and saved to .mat and .csv\n');



%% 4) Plot -every stim per rat
%Individual plots to be run Change as reqd
PL_heatmap_raw=false;
PL_line_trace =true;
PL_mPFC_Regions =true;

% Stim Line conditions
PL_line=true;
PL_regions_line=true;
PL_region_avg =true;

%This runs the plots for each rat
for rat=1:length(T_PRM.RATS)
    rat_no=T_PRM.RATS{rat};
    

    %  Data File Names
    data_filepath=fullfile(T_PRM.STORAGE_DIR,rat_no,sprintf('%s_DATA.mat',rat_no));

    % % Individual LFP Traces
    PRM_Traces_img_dir=fullfile(PRM_ALL_IMG_DIR,'RawTraces_Individual');
    if ~exist(PRM_Traces_img_dir,"dir"); mkdir(PRM_Traces_img_dir); end
    NP2_Plot_LF(rat_no,data_filepath,'PRM_IMG_DIR',PRM_Traces_img_dir,...
       'PL_heatmap_raw',PL_heatmap_raw,...
         'PL_line_trace' ,PL_line_trace,...
         'PL_mPFC_Regions' ,PL_mPFC_Regions);

    % % Responses across all Stim Intensities
     PRM_AllTraces_img_dir=fullfile(PRM_ALL_IMG_DIR,'RawTraces_AcrossStim');
     if ~exist(PRM_AllTraces_img_dir,"dir"); mkdir(PRM_AllTraces_img_dir); end

    % PLot all stim increasing values in one figure
    NP2_Plot_AllStim(rat_no,data_filepath,'PRM_IMG_DIR',PRM_AllTraces_img_dir,...
        'PRM_STORAGE_DIR',T_PRM.STORAGE_DIR,...
        'PL_line',PL_line,'PL_regions_line',PL_regions_line,...
        'PL_region_avg',PL_region_avg);
    
end



%% 5) Get Shank Map (needs manual verification with histology)
LF_Get_Shank_Map(T_PRM);

%% 6) Rat Average Plots 
% Average across LFP of valid stims to get depth traces
NP2_Plot_LFP_DepthTraces(T_PRM);
  
% mPFC Regional Differences (IL and PL)
NP2_Plot_LF_GroupAvg(T_PRM);

% Heatmap Average
NP2_Plot_Heatmap_byDepth(T_PRM);  
