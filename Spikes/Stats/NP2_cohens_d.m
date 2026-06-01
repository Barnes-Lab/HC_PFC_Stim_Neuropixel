function d = NP2_cohens_d(a, b, kind)
% NP2_COHENS_D  Effect size between two samples.
%  kind=paired-mean(a-b) / SD(a-b)
%           =unparied(mean(a) - mean(b)) / pooled SD
% Returns NaN if N < 2 (or unequal N for paired).
% SS 2026

a = a(~isnan(a));
b = b(~isnan(b));

switch lower(kind)
    case 'paired'
        if numel(a) ~= numel(b) || numel(a) < 2; d = NaN; return; end
        diff = a - b;
        d = mean(diff) / std(diff);
    case 'unpaired'
        if numel(a) < 2 || numel(b) < 2; d = NaN; return; end
        s = sqrt(((numel(a)-1)*var(a) + (numel(b)-1)*var(b)) / (numel(a)+numel(b)-2));
        d = (mean(a) - mean(b)) / s;
    otherwise
        error('kind must be ''paired'' or ''unpaired''');
end
end
