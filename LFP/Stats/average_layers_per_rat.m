function data_averaged = average_layers_per_rat(data)
% AVERAGE_LAYERS_PER_RAT - Average L2 and L5 data for each rat in Combined analysis
%
% INPUT:
%   data - data structure with separate L2 and L5 entries per rat
%
% OUTPUT:
%   data_averaged - data structure with one entry per rat (L2 and L5 averaged)

% Get unique rats
unique_rats = unique(data.rat_ids);
n_rats = length(unique_rats);

% Initialize averaged data structure
data_averaged = struct();
data_averaged.amplitude = nan(n_rats, size(data.amplitude, 2), size(data.amplitude, 3), size(data.amplitude, 4));
data_averaged.latency = nan(n_rats, size(data.latency, 2), size(data.latency, 3), size(data.latency, 4));
data_averaged.rat_ids = cell(n_rats, 1);
data_averaged.ages = cell(n_rats, 1);
data_averaged.young_idx = [];
data_averaged.old_idx = [];

% Average L2 and L5 for each rat
for iRat = 1:n_rats
    rat = unique_rats{iRat};
    
    % Find L2 and L5 indices
    L2_idx = find(strcmp(data.rat_ids, rat) & strcmp(data.layers, 'L2'), 1);
    L5_idx = find(strcmp(data.rat_ids, rat) & strcmp(data.layers, 'L5'), 1);
    
    if ~isempty(L2_idx) && ~isempty(L5_idx)
        % Average amplitude and latency across L2 and L5
        data_averaged.amplitude(iRat, :, :, :) = nanmean(cat(1, ...
            reshape(data.amplitude(L2_idx, :, :, :), [1 size(data.amplitude, 2:4)]), ...
            reshape(data.amplitude(L5_idx, :, :, :), [1 size(data.amplitude, 2:4)])), 1);
        
        data_averaged.latency(iRat, :, :, :) = nanmean(cat(1, ...
            reshape(data.latency(L2_idx, :, :, :), [1 size(data.latency, 2:4)]), ...
            reshape(data.latency(L5_idx, :, :, :), [1 size(data.latency, 2:4)])), 1);
        
        % Store rat info
        data_averaged.rat_ids{iRat} = rat;
        data_averaged.ages{iRat} = data.ages{L2_idx};
        
        % Build age indices
        if strcmp(data.ages{L2_idx}, 'Y')
            data_averaged.young_idx(end+1) = iRat;
        else
            data_averaged.old_idx(end+1) = iRat;
        end
    end
end

end
