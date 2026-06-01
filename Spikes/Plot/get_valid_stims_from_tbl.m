function valid = get_valid_stims_from_tbl(tbl, rat_id, sess_name, sorted_current_uA, varargin)
% GET_VALID_STIMS_FROM_TBL  Return indices of stims with Keep != 'N' for one
% (rat, session). Defaults to 'Keep' (capital K); auto-falls back to 'keep'
% if 'Keep' is not in the table. Replaces the duplicated CSV-lookup logic
% across the NP2 plotting and analysis functions.
%
% Example:
%   valid = get_valid_stims_from_tbl(tbl, '10987', 'iHC', STIM.sorted_current_uA);
%   valid = get_valid_stims_from_tbl(tbl, ..., 'PRM_keep_col', 'keep');
%
% INPUT
%   tbl                table  ALL_STIM_TIMES_KEEP_CHECK.csv rows
%   rat_id             char | string
%   sess_name          char | string  e.g. 'iHC' | 'vHC'
%   sorted_current_uA  [1xN]  stim currents in the same order as STIM
%   varargin           name-value:
%     'PRM_keep_col'   char   force a specific column ('Keep' | 'keep'),
%                             otherwise auto-detected
%
% OUTPUT
%   valid              [1xK]  indices into sorted_current_uA to keep
%
% SS 2026

PRM_keep_col = '';
Extract_varargin;

% Resolve column: explicit override > 'Keep' > 'keep'
if isempty(PRM_keep_col)
    if ismember('Keep', tbl.Properties.VariableNames)
        keep_col = 'Keep';
    elseif ismember('keep', tbl.Properties.VariableNames)
        keep_col = 'keep';
    else
        error('Neither ''Keep'' nor ''keep'' column found in stim table.');
    end
else
    keep_col = PRM_keep_col;
end

rows = tbl(strcmp(string(tbl.Rat), string(rat_id)) & ...
           strcmp(tbl.Session,    string(sess_name)), :);

ex = false(1, numel(sorted_current_uA));
for k = 1:numel(sorted_current_uA)
    idx = find(rows.Current_uA == round(sorted_current_uA(k)), 1);
    if ~isempty(idx) && strcmp(rows.(keep_col){idx}, 'N')
        ex(k) = true;
    end
end
valid = find(~ex);
end
