function sr = empty_shank_result()
% Default empty shank-result struct used to pre-allocate 
%
% SS 2026

sr = struct();
sr.shank_idx         = NaN;
sr.layer             = '';
sr.sess_name         = '';
sr.channel_mode      = '';
sr.stim_currents_all = [];
sr.valid_stims       = [];
sr.x                 = [];
sr.y                 = [];
sr.has_data          = false;
sr.is_flat           = false;
sr.fit               = struct('I80',NaN,'Amin',NaN,'Amax',NaN,'I50',NaN, ...
                              'k',NaN,'x0',NaN,'R2',NaN,'Method','none');
sr.Status            = 'no_data';

end
