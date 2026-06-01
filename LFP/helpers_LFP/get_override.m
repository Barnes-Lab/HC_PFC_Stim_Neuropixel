function ovr = get_override(OVR_struct, rat_field, sess_name)

% Looks up T_PRM.OVR.RAT.SESS and returns the override
% struct -Used by RUN_LF_I80_DETECTION 



ovr = [];
if ~isstruct(OVR_struct);                       return; end
if ~isfield(OVR_struct, rat_field);             return; end
sub = OVR_struct.(rat_field);
if ~isfield(sub, sess_name);                    return; end
ovr = sub.(sess_name);
end
