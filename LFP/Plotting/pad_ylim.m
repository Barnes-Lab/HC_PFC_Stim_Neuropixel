function ylim_vec = pad_ylim(vals, pad_low, pad_high)
% Compute padded y-limits from a value array.
%
%
% SS 2026

if nargin < 2 || isempty(pad_low);  pad_low  = 0.15; end
if nargin < 3 || isempty(pad_high); pad_high = 0.25; end

v = vals(isfinite(vals));
if isempty(v); ylim_vec = [0, 1]; return; end

vmin = min(v); vmax = max(v);
if vmin == vmax
    if vmin == 0; ylim_vec = [-0.1, 0.1];
    else;         ylim_vec = [vmin - abs(vmin)*0.1, vmax + abs(vmax)*0.1];
    end
    return;
end
r        = vmax - vmin;
ylim_vec = [vmin - pad_low*r, vmax + pad_high*r];
end
