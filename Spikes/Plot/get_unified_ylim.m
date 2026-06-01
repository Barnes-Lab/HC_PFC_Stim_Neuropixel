function ylim_vals = get_unified_ylim(data_cell, padding_factor)
% Unified y-axis limits across multiple datasets.



if nargin < 2; padding_factor = 0.1; end

all_data = [];
for i = 1:length(data_cell)
    if iscell(data_cell{i})
        for j = 1:length(data_cell{i})
            if ~isempty(data_cell{i}{j}); all_data = [all_data; data_cell{i}{j}(:)]; end
        end
    elseif ~isempty(data_cell{i})
        all_data = [all_data; data_cell{i}(:)];
    end
end
all_data = all_data(isfinite(all_data));

if isempty(all_data)
    ylim_vals = [0, 1];
    warning('No valid data for ylim calculation, using default [0, 1]');
    return;
end

ymin = min(all_data);
ymax = max(all_data);
if ymin == ymax
    if ymin == 0; ylim_vals = [-0.1, 0.1];
    else;         ylim_vals = [ymin - abs(ymin)*0.1, ymax + abs(ymax)*0.1];
    end
    return;
end

range_val = ymax - ymin;
ylim_vals = [ymin - padding_factor*range_val, ymax + padding_factor*range_val];
end
