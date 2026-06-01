function result = try_fit_cascade(x, y, maybe_mask, varargin)
% Cascade fit: tries logistic, logistic-Amin-zero, exponential saturation
%
% Stops at the FIRST method that satisfies:
%   fit_success && R2 >= R2_thr && I80 in PRM_I80_range
%
% DEPENDS  fit_logistic_constrained.m, fit_saturation_exp.m, data_driven_i80.m
%
% SS 2026

%% PARAMETERS TO CHANGE
PRM_R2_min       = 0.50;
PRM_R2_forced    = 0.20;
PRM_I80_range    = [50, 900];
PRM_force_method = '';
Extract_varargin;

%% Force-method shortcut
if ~isempty(PRM_force_method)
    result = run_force_method(x, y, PRM_force_method);
    return;
end

%% All stims
result = try_all_methods(x, y, '', PRM_R2_min, PRM_I80_range);
if passed(result); return; end


if any(maybe_mask)
    xm = x(~maybe_mask); ym = y(~maybe_mask);
    if numel(xm) >= 4
        r2 = try_all_methods(xm, ym, '_no_maybe', PRM_R2_min, PRM_I80_range);
        if passed(r2); result = r2; return; end
    end
end

%
result = try_all_methods(x, y, '_forced', PRM_R2_forced, PRM_I80_range);
if passed(result); return; end

%
if any(maybe_mask)
    xm = x(~maybe_mask); ym = y(~maybe_mask);
    if numel(xm) >= 4
        r2 = try_all_methods(xm, ym, '_forced_no_maybe', PRM_R2_forced, PRM_I80_range);
        if passed(r2); result = r2; return; end
    end
end


result        = default_result();
[I80_dd, plateau_dd] = data_driven_i80(x, y);
result.I80    = I80_dd;
result.Amax   = plateau_dd;
result.Amin   = 0;
result.Method = 'data_driven';

end

% Inline functions short
function r = try_all_methods(x, y, suffix, R2_thr, I80_rng)
r = default_result();

f = fit_logistic_constrained(x, y);
if check_fit(f, R2_thr, I80_rng)
    r = pack_logistic(f, ['logistic' suffix]); return;
end

f = fit_logistic_constrained(x, y, 'PRM_force_amin_zero', true);
if check_fit(f, R2_thr, I80_rng)
    r = pack_logistic(f, ['logistic_A0' suffix]); return;
end

f = fit_saturation_exp(x, y);
if check_fit(f, R2_thr, I80_rng)
    r = pack_expsat(f, ['exp_sat' suffix]); return;
end
end

function ok = check_fit(f, R2_thr, I80_rng)
ok = isstruct(f) && f.fit_success && ~isnan(f.R2) && f.R2 >= R2_thr && ...
     ~isnan(f.I80) && f.I80 >= I80_rng(1) && f.I80 <= I80_rng(2);
end

function r = default_result()
r = struct('I80', NaN, 'Amin', NaN, 'Amax', NaN, 'I50', NaN, ...
    'k', NaN, 'x0', NaN, 'R2', NaN, 'Method', 'none');
end

function r = pack_logistic(f, method)
r        = default_result();
r.I80    = f.I80;  r.Amin = f.Amin;  r.Amax = f.Amax;
r.I50    = f.I50;  r.k    = f.k;     r.R2   = f.R2;
r.Method = method;
end

function r = pack_expsat(f, method)
r        = default_result();
r.I80    = f.I80;  r.Amin = f.Amin;  r.Amax = f.Amax;
r.I50    = f.I50;  r.k    = f.k;     r.x0   = f.x0;  r.R2 = f.R2;
r.Method = method;
end

function tf = passed(result)
tf = ~strcmp(result.Method, 'none') && ~strcmp(result.Method, 'data_driven');
end
