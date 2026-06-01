function rat_means = calculate_rat_means(rat_ids, values)
% CALCULATE_RAT_MEANS  Mean value per unique rat. 
% Rat id can be numeric or cell (with char)
%
% %
% SS 2026

if iscell(rat_ids)
    rat_ids_numeric = zeros(length(rat_ids), 1);
    for i = 1:length(rat_ids)
        element = rat_ids{i};
        while iscell(element)
            if isempty(element); element = NaN; break; end
            element = element{1};
        end
        if isnumeric(element);                   rat_ids_numeric(i) = element;
        elseif ischar(element) || isstring(element); rat_ids_numeric(i) = str2double(element);
        else;                                    rat_ids_numeric(i) = NaN;
        end
    end
else
    rat_ids_numeric = rat_ids(:);
end

unique_rats = unique(rat_ids_numeric(~isnan(rat_ids_numeric)));
rat_means   = zeros(length(unique_rats), 1);
for r = 1:length(unique_rats)
    rat_means(r) = nanmean(values(rat_ids_numeric == unique_rats(r)));
end
end
