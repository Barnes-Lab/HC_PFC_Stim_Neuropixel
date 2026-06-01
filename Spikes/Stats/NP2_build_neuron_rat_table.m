function T = NP2_build_neuron_rat_table(data, neuron_type, sess_types, recruit_thresh)
% NP2_BUILD_NEURON_RAT_TABLE  Aggregate per-neuron data to per-(rat, sess,
% region, window) rows. 
% NEURAL METRICS HERE delta_FR (mean trial-avg dFR, Hz),
% recruitment (% neurons with mean dFR > thresh), 
% reliability_per_trial (mean per-neuron % trials > thresh).
%
%
% INPUT
%   data-output of collect_neuron_data_with_I80
%   neuron_type-'INT' | 'PYR'
%   sess_types {'iHC','vHC'} 
%   recruit_thresh Hz 
%
% OUTPUT
%   table: Rat, Age, Sess, Region, Window, delta_FR, recruitment,
%             reliability_per_trial, N_neurons
%
% SS 2026

if nargin < 4 || isempty(recruit_thresh); recruit_thresh = 1; end

regions = {'PL', 'IL'};
windows = {'mono', 'poly'};
ages    = {'Y', 'O'};

rows = {};
for s = 1:numel(sess_types)
    sess = sess_types{s};
    for ai = 1:2
        a = ages{ai};
        d = data.(sess).(a).(neuron_type);
        if isempty(d.rat); continue; end

        rats_u = unique(d.rat);
        for ri = 1:numel(rats_u)
            rat      = rats_u{ri};
            rat_mask = strcmp(d.rat, rat);
            for rg = 1:2
                region    = regions{rg};
                cell_mask = rat_mask & strcmp(d.region, region);
                n         = sum(cell_mask);
                if n == 0; continue; end

                for w = 1:2
                    if w == 1
                        vals = d.win1(cell_mask); rel = d.win1_reliability(cell_mask);
                    else
                        vals = d.win2(cell_mask); rel = d.win2_reliability(cell_mask);
                    end
                    rows(end+1, :) = {rat, a, sess, region, windows{w}, ...
                        mean(vals, 'omitnan'), ...
                        100 * sum(vals > recruit_thresh) / n, ...
                        mean(rel, 'omitnan'), n};
                end
            end
        end
    end
end

T = cell2table(rows, 'VariableNames', {'Rat','Age','Sess','Region','Window', ...
    'delta_FR','recruitment','reliability_per_trial','N_neurons'});
T.Age    = categorical(T.Age,    {'Y','O'});
T.Sess   = categorical(T.Sess,   sess_types);
T.Region = categorical(T.Region, {'PL','IL'});
T.Window = categorical(T.Window, {'mono','poly'});
T.Rat    = categorical(T.Rat);
end
