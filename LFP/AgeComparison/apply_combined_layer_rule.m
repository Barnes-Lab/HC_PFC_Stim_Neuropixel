function data_pool = apply_combined_layer_rule(data)
% Collapse combined (L2+L5) data to one row per rat
%   PL: L5 only -
%   IL : mean of L2 and L5 per rat 
% SS 2026

unique_rats = unique(data.rat_ids, 'stable');
n_rats      = numel(unique_rats);
sz          = size(data.amplitude);

amp_pool = nan([n_rats, sz(2:end)]);
lat_pool = nan([n_rats, sz(2:end)]);
ages     = cell(n_rats, 1);

% Cant vecotrize
for iRat = 1:n_rats
    rat    = unique_rats{iRat};
    L2_idx = find(strcmp(data.rat_ids, rat) & strcmp(data.layers, 'L2'), 1);
    L5_idx = find(strcmp(data.rat_ids, rat) & strcmp(data.layers, 'L5'), 1);

    if ~isempty(L5_idx);     ages{iRat} = data.ages{L5_idx};
    elseif ~isempty(L2_idx); ages{iRat} = data.ages{L2_idx};
    end

    % L5 only for PL
    use_L2_PL = false;
    if isempty(L5_idx)
        L5_missing_PL = true;
    else
        L5_missing_PL = any(isnan(reshape(data.latency(L5_idx, :, 1, :), [], 1)));
    end
    if L5_missing_PL && ~isempty(L2_idx)
        L2_has_PL = any(~isnan(reshape(data.latency(L2_idx, :, 1, :), [], 1)));
        if L2_has_PL
            use_L2_PL = true;
            fprintf('  [INFO] Rat %s: L5 PL incomplete -> using L2 PL\n', rat);
        end
    end

    if use_L2_PL
        amp_pool(iRat, :, 1, :) = data.amplitude(L2_idx, :, 1, :);
        lat_pool(iRat, :, 1, :) = data.latency  (L2_idx, :, 1, :);
    elseif ~isempty(L5_idx)
        amp_pool(iRat, :, 1, :) = data.amplitude(L5_idx, :, 1, :);
        lat_pool(iRat, :, 1, :) = data.latency  (L5_idx, :, 1, :);
    end

    % IL boht
    if ~isempty(L2_idx) && ~isempty(L5_idx)
        amp_pool(iRat, :, 2, :) = nanmean(cat(1, ...
            data.amplitude(L2_idx, :, 2, :), data.amplitude(L5_idx, :, 2, :)), 1);
        lat_pool(iRat, :, 2, :) = nanmean(cat(1, ...
            data.latency  (L2_idx, :, 2, :), data.latency  (L5_idx, :, 2, :)), 1);
    elseif ~isempty(L5_idx)
        amp_pool(iRat, :, 2, :) = data.amplitude(L5_idx, :, 2, :);
        lat_pool(iRat, :, 2, :) = data.latency  (L5_idx, :, 2, :);
    elseif ~isempty(L2_idx)
        amp_pool(iRat, :, 2, :) = data.amplitude(L2_idx, :, 2, :);
        lat_pool(iRat, :, 2, :) = data.latency  (L2_idx, :, 2, :);
    end
end

data_pool.amplitude = amp_pool;
data_pool.latency   = lat_pool;
data_pool.rat_ids   = unique_rats(:);
data_pool.ages      = ages;
data_pool.young_idx = find(strcmp(ages, 'Y'));
data_pool.old_idx   = find(strcmp(ages, 'O'));

end
