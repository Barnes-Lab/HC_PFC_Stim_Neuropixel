function T_filt = filter_stim_times_to_I80(ALL_STIM_TIMES, I80_DATA, sess_types)
% FILTER_STIM_TIMES_TO_I80  Copy of ALL_STIM_TIMES containing only the 5
% I80-selected stim rows per (rat, sess). Same columns as input.
% Drop-in replacement for the full table when calling PSTH functions in
% I80 mode.
% SS 2026

cur_col = '';
for cc = {'Current_uA', 'Intensity_uA', 'Stim_uA'}
    if ismember(cc{1}, ALL_STIM_TIMES.Properties.VariableNames); cur_col = cc{1}; break; end
end
if isempty(cur_col); error('Stim-current column not found in ALL_STIM_TIMES'); end

rat_col   = string(ALL_STIM_TIMES.Rat);
rats_uniq = unique(rat_col);
keep      = false(height(ALL_STIM_TIMES), 1);

for ri = 1:numel(rats_uniq)
    rat_id = char(rats_uniq(ri));
    for si = 1:numel(sess_types)
        sess    = sess_types{si};
        rs_mask = (rat_col == rats_uniq(ri)) & strcmp(ALL_STIM_TIMES.Session, sess);
        if ~any(rs_mask); continue; end

        T_sess = ALL_STIM_TIMES(rs_mask, :);
        [ev_times, ~, ~] = select_I80_event_times_neuron(I80_DATA, rat_id, sess, T_sess, cur_col);
        if isempty(ev_times); continue; end

        sub_idx = find(rs_mask);
        for ei = 1:numel(ev_times)
            local = abs(T_sess.Time_s - ev_times(ei)) < 1e-6;
            keep(sub_idx(local)) = true;
        end
    end
end
T_filt = ALL_STIM_TIMES(keep, :);
end
