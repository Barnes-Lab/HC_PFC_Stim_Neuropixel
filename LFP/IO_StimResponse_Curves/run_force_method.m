function r = run_force_method(x, y, method)
% No R2 / I80 range caps — caller decides what to do with output
r = default_result();
switch lower(method)
    case 'auto'
        % Try logistic, logistic_A0, exp_sat — return highest-R2 fit
        candidates = {};
        f = fit_logistic_constrained(x, y);
        if isstruct(f) && f.fit_success
            candidates{end+1} = pack_logistic(f, 'logistic_auto'); %#ok<AGROW>
        end
        f = fit_logistic_constrained(x, y, 'PRM_force_amin_zero', true);
        if isstruct(f) && f.fit_success
            candidates{end+1} = pack_logistic(f, 'logistic_A0_auto'); %#ok<AGROW>
        end
        f = fit_saturation_exp(x, y);
        if isstruct(f) && f.fit_success
            candidates{end+1} = pack_expsat(f, 'exp_sat_auto'); %#ok<AGROW>
        end
        if isempty(candidates); return; end
        R2s = cellfun(@(c) c.R2, candidates);
        [~, best] = max(R2s);
        r = candidates{best};
        r.Method = [r.Method '_forced_method'];

    case 'logistic'
        f = fit_logistic_constrained(x, y);
        if isstruct(f) && f.fit_success
            r = pack_logistic(f, 'logistic_forced_method');
        end
    case 'logistic_a0'
        f = fit_logistic_constrained(x, y, 'PRM_force_amin_zero', true);
        if isstruct(f) && f.fit_success
            r = pack_logistic(f, 'logistic_A0_forced_method');
        end
    case 'exp_sat'
        f = fit_saturation_exp(x, y);
        if isstruct(f) && f.fit_success
            r = pack_expsat(f, 'exp_sat_forced_method');
        end
    case 'data_driven'
        [I80_dd, plat]  = data_driven_i80(x, y);
        r.I80           = I80_dd;
        r.Amax          = plat;
        r.Amin          = 0;
        r.Method        = 'data_driven_forced';
    otherwise
        error('try_fit_cascade: unknown force method "%s"', method);
end
end
