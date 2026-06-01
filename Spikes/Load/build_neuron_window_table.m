function W = build_neuron_window_table(rats, sess_types, regions, mode, defaults, ...
    manual_csv, fepsp_csv, fepsp_win_type)
% BUILD_NEURON_WINDOW_TABLE  Per-(rat, sess, region) mono/poly window table
% for spike analysis.
% 
% SS 2026

n   = numel(rats) * numel(sess_types) * numel(regions);
Rat = cell(n, 1); Sess = cell(n, 1); Region = cell(n, 1);
ms  = zeros(n, 1); me = zeros(n, 1); ps = zeros(n, 1); pe = zeros(n, 1);

i = 0;
for ri = 1:numel(rats)
    for si = 1:numel(sess_types)
        for rg = 1:numel(regions)
            i = i + 1;
            Rat{i}    = rats{ri};
            Sess{i}   = sess_types{si};
            Region{i} = regions{rg};
            ms(i) = defaults.mono(1); me(i) = defaults.mono(2);
            ps(i) = defaults.poly(1); pe(i) = defaults.poly(2);
        end
    end
end
W = table(Rat, Sess, Region, ms, me, ps, pe, ...
    'VariableNames', {'Rat','Sess','Region','Mono_Start','Mono_End','Poly_Start','Poly_End'});

switch mode
    case 'default';    return;
    case 'manual';     W = override_manual(W, manual_csv);
    case 'fepsp_peak'; W = override_fepsp(W, fepsp_csv, fepsp_win_type);
    otherwise;         error('Unknown mode: %s', mode);
end
end
