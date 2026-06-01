function [data, stim_report] = collect_neuron_data_with_I80( ...
    all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, ...
    window_table, bsl_ms, PL_depth, sess_types, neuron_types, ...
    I80_DATA, recruit_thresh)
% Per-neuron responses locked to 5 I80-selectedstims (simlar to LFP). 
% Windows looked up per (rat, sess, region) from
% window_table. Three response windows per neuron:
%   win1- mono-[Mono_Start, Mono_End]
%   win2 - poly -[Poly_Start, Poly_End]
%   win_combined -[Mono_Start, Poly_End]
-
%
% INPUT
%   all_SP_filtered, ALL_STIM_TIMES, neuron_classifications, I80_DATA
%    build_neuron_window_table
%   bsl_ms 
%   PL_depth  default is 3600uM
%   sess_types- {'iHC','vHC'}
%   neuron_types  - {'INT','PYR'}
%   recruit_thresh  - Hz
%
% 
% SS 2026

if nargin < 10 || isempty(recruit_thresh); recruit_thresh = 1; end
warned_degenerate = false;

%% Init data struct
data = struct();
for s = 1:numel(sess_types)
    sess = sess_types{s};
    for age = {'Y','O'}
        a = age{1};
        for t = 1:numel(neuron_types)
            ty = neuron_types{t};
            data.(sess).(a).(ty) = struct( ...
                'rat', {{}}, 'region', {{}}, ...
                'win1', [], 'win2', [], 'win_combined', [], ...
                'win1_reliability', [], 'win2_reliability', [], ...
                'win1_5of5', logical([]), 'win2_5of5', logical([]));
        end
    end
end
report_entries = {};

% Loop neurons
for u = 1:numel(all_SP_filtered)
    SP = all_SP_filtered(u);
    nt = neuron_classifications.type{u};
    if strcmp(nt, 'UNCLASSIFIED'); continue; end

    rat = SP.rat;
    age = SP.age;
    if isfield(SP, 'neuropixels_depth_uM'); dep = SP.neuropixels_depth_uM;
    elseif isfield(SP, 'depth_uM');         dep = SP.depth_uM;
    else; continue;
    end
    region = 'IL'; if dep < PL_depth; region = 'PL'; end

    spk_s = double(SP.t_uS) / 1e6;

    for s = 1:numel(sess_types)
        sess = sess_types{s};
        [win1_ms, win2_ms] = lookup_neuron_windows(window_table, rat, sess, region);
        if isempty(win1_ms); continue; end

        win_c_ms = [win1_ms(1), win2_ms(2)];
        if win_c_ms(2) <= win_c_ms(1)
            if ~warned_degenerate
                warning('win_combined degenerate (mono_start >= poly_end) for %s %s %s -- setting NaN', ...
                    rat, sess, region);
                warned_degenerate = true;
            end
            win_c_ms = [];
        end

        T = ALL_STIM_TIMES(strcmp(string(ALL_STIM_TIMES.Rat), rat) & ...
                           strcmp(ALL_STIM_TIMES.Session, sess), :);
        if isempty(T); continue; end

        cur_col = '';
        for cc = {'Current_uA','Intensity_uA','Stim_uA'}
            if ismember(cc{1}, T.Properties.VariableNames); cur_col = cc{1}; break; end
        end
        if isempty(cur_col); continue; end

        [event_times, currents_5, method] = select_I80_event_times_neuron( ...
            I80_DATA, rat, sess, T, cur_col);
        if isempty(event_times); continue; end

        resp = compute_response_metrics(spk_s, event_times, ...
            win1_ms, win2_ms, win_c_ms, bsl_ms, recruit_thresh);

        d = data.(sess).(age).(nt);
        d.rat{end+1, 1}              = rat;
        d.region{end+1, 1}           = region;
        d.win1(end+1, 1)             = resp.win1;
        d.win2(end+1, 1)             = resp.win2;
        d.win_combined(end+1, 1)     = resp.win_combined;
        d.win1_reliability(end+1, 1) = resp.win1_reliability;
        d.win2_reliability(end+1, 1) = resp.win2_reliability;
        d.win1_5of5(end+1, 1)        = resp.win1_5of5;
        d.win2_5of5(end+1, 1)        = resp.win2_5of5;
        data.(sess).(age).(nt) = d;

        report_entries(end+1, :) = {rat, age, sess, region, nt, ...
            sprintf('[%g %g]/[%g %g]', win1_ms, win2_ms), ...
            numel(event_times), method, ...
            sprintf('%.0f-%.0f', min(currents_5), max(currents_5))};
    end
end

stim_report = cell2table(report_entries, 'VariableNames', ...
    {'Rat','Age','Session','Region','Type','Windows_ms','N_Stims','Method','Range_uA'});
end


function resp = compute_response_metrics(spk_s, ev_s, win1_ms, win2_ms, win_c_ms, bsl_ms, thresh)
% Per-trial dFR for mono, poly, combined windows.
w1 = win1_ms / 1000; w2 = win2_ms / 1000; b = bsl_ms / 1000;
dw1 = diff(w1); dw2 = diff(w2); db = diff(b);
n   = numel(ev_s);

t1 = zeros(n, 1); t2 = zeros(n, 1); tc = nan(n, 1);
if ~isempty(win_c_ms); wc = win_c_ms / 1000; dwc = diff(wc);
else;                  wc = [];               dwc = NaN;
end

for e = 1:n
    t0   = ev_s(e);
    bsl  = sum(spk_s >= t0 + b(1)  & spk_s < t0 + b(2))  / db;
    r1   = sum(spk_s >= t0 + w1(1) & spk_s < t0 + w1(2)) / dw1;
    r2   = sum(spk_s >= t0 + w2(1) & spk_s < t0 + w2(2)) / dw2;
    t1(e) = r1 - bsl;
    t2(e) = r2 - bsl;
    if ~isempty(wc)
        rc    = sum(spk_s >= t0 + wc(1) & spk_s < t0 + wc(2)) / dwc;
        tc(e) = rc - bsl;
    end
end

resp.win1             = mean(t1);
resp.win2             = mean(t2);
resp.win_combined     = mean(tc, 'omitnan');
resp.win1_reliability = 100 * sum(t1 > thresh) / n;
resp.win2_reliability = 100 * sum(t2 > thresh) / n;
resp.win1_5of5        = all(t1 > thresh);
resp.win2_5of5        = all(t2 > thresh);
end
