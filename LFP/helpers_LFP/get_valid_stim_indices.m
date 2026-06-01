function valid_stims = get_valid_stim_indices(allStimTbl, rat_id, sess_name, ...
                                                sorted_current_uA)

% Returns valid stim indices
%
% SS 2026

nStims      = numel(sorted_current_uA);
valid_stims = 1:nStims;

if isempty(allStimTbl); return; end

rows = allStimTbl( ...
    strcmp(string(allStimTbl.Rat),     string(rat_id)) & ...
    strcmp(string(allStimTbl.Session), string(sess_name)), :);
if isempty(rows); return; end

varNames = rows.Properties.VariableNames;
keep_var = '';
if     any(strcmp(varNames, 'Keep')); keep_var = 'Keep';
elseif any(strcmp(varNames, 'keep')); keep_var = 'keep';
end
if isempty(keep_var); return; end

keep_col = string(rows.(keep_var));
n_use    = min(numel(keep_col), nStims);

excluded = false(1, nStims);
for k = 1:n_use
    if strcmpi(keep_col(k), 'N'); excluded(k) = true; end
end
valid_stims = find(~excluded);

end
