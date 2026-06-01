function sess_name = peek_session_name(rat_struct, sess_idx)
% Returns the session name ('iHC' or 'vHC') for the given session
% SS 2026

sess_name = '';
flds = fieldnames(rat_struct);
flds = flds(startsWith(flds, 'Shank'));
for i = 1:numel(flds)
    sh = rat_struct.(flds{i});
    if ~isfield(sh, 'sessions') || numel(sh.sessions) < sess_idx; continue; end
    if isfield(sh.sessions(sess_idx), 'name') && ~isempty(sh.sessions(sess_idx).name)
        sess_name = sh.sessions(sess_idx).name; return;
    end
end
end
