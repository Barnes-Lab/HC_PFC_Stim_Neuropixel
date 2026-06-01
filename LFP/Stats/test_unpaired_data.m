function [p, h] = test_unpaired_data(data1, data2)
% TEST_UNPAIRED_DATA - Mann-Whitney U test for unpaired data
%
% Returns p-value and hypothesis test result
% Requires at least 3 valid observations in each group

valid1 = ~isnan(data1);
valid2 = ~isnan(data2);

if sum(valid1) >= 3 && sum(valid2) >= 3
    [p, h] = ranksum(data1(valid1), data2(valid2));
else
    p = NaN;
    h = 0;
end

end
