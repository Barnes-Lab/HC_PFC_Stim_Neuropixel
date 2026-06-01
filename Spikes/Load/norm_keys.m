function [rat, sess, region] = norm_keys(rat_col, sess_col, region_col)
% NORM_KEYS  Normalize rat/sess/region table columns to char or string



if isnumeric(rat_col); rat = arrayfun(@num2str, rat_col, 'UniformOutput', false);
else;                  rat = strtrim(string(rat_col));
end
sess   = strtrim(string(sess_col));
region = strtrim(string(region_col));
end
