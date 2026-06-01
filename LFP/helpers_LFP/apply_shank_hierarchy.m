function decision = apply_shank_hierarchy(shank_results, L5_shank, override, varargin)
%
% Decides which shank's I80 to use for the (rat, session) compares them for
% diff fit. Not used in the end
%
%
% SS 2026

%% PARAMETERS
PRM_I80_diff_thresh_uA = 100;
PRM_flat_force_uA      = 200;
Extract_varargin;

decision = struct( ...
    'chosen_shank', NaN, 'chosen_source', 'no_usable_shank', ...
    'I80', NaN, 'L5_I80', NaN, ...
    'Best_R2_Shank', NaN, 'L5_vs_Best_diff_uA', NaN);

n        = numel(shank_results);
ok_class = {'ok','ok_no_maybe','ok_forced','ok_forced_no_maybe','ok_force_method'};
is_ok    = false(n,1);
R2_arr   = nan(n,1);
I80_arr  = nan(n,1);

for i = 1:n
    sr = shank_results(i);
    if any(strcmp(sr.Status, ok_class))
        is_ok(i)   = true;
        R2_arr(i)  = sr.fit.R2;
        I80_arr(i) = sr.fit.I80;
    end
end

[best_R2, best_idx] = max(R2_arr, [], 'omitnan');
if isnan(best_R2); best_idx = NaN; end
decision.Best_R2_Shank = best_idx;

L5_valid = ~isnan(L5_shank) && L5_shank >= 1 && L5_shank <= n;
L5_ok    = L5_valid && is_ok(L5_shank);
if L5_ok; decision.L5_I80 = I80_arr(L5_shank); end

% Manual (always used)
if ~isempty(override) && isstruct(override) && isfield(override, 'shank')
    sh = override.shank;
    if sh >= 1 && sh <= n
        sr = shank_results(sh);
        decision.chosen_shank = sh;

        if isfield(override, 'force_flat') && override.force_flat
            decision.chosen_source = 'manual_override_flat';
            decision.I80           = PRM_flat_force_uA;
        else
            decision.chosen_source = 'manual_override';
            decision.I80           = sr.fit.I80;
            % If override fit returned NaN, fall back to data-driven on its data
            if isnan(decision.I80) && sr.has_data
                [I80_dd, ~]  = data_driven_i80(sr.x, sr.y);
                decision.I80 = I80_dd;
            end
        end
        return;
    end
end

% % L5 ok
% if L5_ok
%     if ~isnan(best_idx) && best_idx ~= L5_shank
%         decision.L5_vs_Best_diff_uA = abs(I80_arr(L5_shank) - I80_arr(best_idx));
%         if decision.L5_vs_Best_diff_uA > PRM_I80_diff_thresh_uA
%             decision.chosen_shank  = best_idx;
%             decision.chosen_source = 'best_R2_diverged';
%             decision.I80           = I80_arr(best_idx);
%         else
%             decision.chosen_shank  = L5_shank;
%             decision.chosen_source = 'user_L5';
%             decision.I80           = I80_arr(L5_shank);
%         end
%     else
%         decision.chosen_shank       = L5_shank;
%         decision.chosen_source      = 'user_L5';
%         decision.I80                = I80_arr(L5_shank);
%         decision.L5_vs_Best_diff_uA = 0;
%     end
%     return;
% end
% 
% %
% if ~isnan(best_idx)
%     decision.chosen_shank  = best_idx;
%     decision.chosen_source = 'best_R2_L5_failed';
%     decision.I80           = I80_arr(best_idx);
%     return;
% end
% 
% %
% if L5_valid
%     sr = shank_results(L5_shank);
%     if sr.is_flat
%         decision.chosen_shank  = L5_shank;
%         decision.chosen_source = 'forced_flat';
%         decision.I80           = PRM_flat_force_uA;
%         return;
%     elseif sr.has_data && ~isnan(sr.fit.I80)
%         decision.chosen_shank  = L5_shank;
%         decision.chosen_source = 'data_driven';
%         decision.I80           = sr.fit.I80;
%         return;
%     end
% end
% 
% for i = 1:n
%     if shank_results(i).is_flat
%         decision.chosen_shank  = i;
%         decision.chosen_source = 'forced_flat_other_shank';
%         decision.I80           = PRM_flat_force_uA;
%         return;
%     end
% end
% 
% for i = 1:n
%     sr = shank_results(i);
%     if sr.has_data && ~isnan(sr.fit.I80)
%         decision.chosen_shank  = i;
%         decision.chosen_source = 'data_driven_other_shank';
%         decision.I80           = sr.fit.I80;
%         return;
%     end
% end

end
