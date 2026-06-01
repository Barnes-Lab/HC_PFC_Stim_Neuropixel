function [L2_shank, L5_shank] = lookup_layer_shanks(rat_id, shankMap)
% Returns L2 and L5 shank numbers from ShankMap.csv for a given rat.
% Both NaN if shankMap is empty or rat is not present in the map.
%
% SS 2026

L2_shank = NaN;
L5_shank = NaN;

if isempty(shankMap); return; end

mi = find(shankMap.Rat == str2double(rat_id), 1);
if isempty(mi); return; end

L2_shank = shankMap.L2(mi);
L5_shank = shankMap.L5(mi);

end
