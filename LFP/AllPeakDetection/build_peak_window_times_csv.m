function build_peak_window_times_csv(SR_DATA, T_VEC, T_PRM, save_path)
% BUILD_PEAK_WINDOW_TIMES_CSV  Per-rat peak window times into CSV.
%

% OUTPUT:
%   CSV with columns:
%     Rat, Age, Session, Region, Window,
%     Peak_Latency_L2_ms, Peak_Amplitude_L2_mV,
%     Peak_Latency_L5_ms, Peak_Amplitude_L5_mV,
%     HalfMax_Start_ms, HalfMax_End_ms, HalfMax_Duration_ms,
%     Zero_Start_ms, Zero_End_ms, Zero_Duration_ms
%
% SS 2026

regions  = {'PL', 'IL'};
windows  = {'mono', 'poly'};
win_lbls = {'Mono', 'Poly'};
sessions = {'iHC', 'vHC'};

n_rats   = numel(T_PRM.RATS);
n_rows   = n_rats * numel(sessions) * numel(regions) * numel(windows);

% Preallocate
Rat                    = strings(n_rows, 1);
Age                    = strings(n_rows, 1);
Session                = strings(n_rows, 1);
Region                 = strings(n_rows, 1);
Window                 = strings(n_rows, 1);
Peak_Latency_L2_ms     = nan(n_rows, 1);
Peak_Amplitude_L2_mV   = nan(n_rows, 1);
Peak_Latency_L5_ms     = nan(n_rows, 1);
Peak_Amplitude_L5_mV   = nan(n_rows, 1);
HalfMax_Start_ms       = nan(n_rows, 1);
HalfMax_End_ms         = nan(n_rows, 1);
HalfMax_Duration_ms    = nan(n_rows, 1);
Zero_Start_ms          = nan(n_rows, 1);
Zero_End_ms            = nan(n_rows, 1);
Zero_Duration_ms       = nan(n_rows, 1);

row = 0;
for r = 1:n_rats
    rat_id = T_PRM.RATS{r};
    age    = T_PRM.AGES{r};
    ratFld = sprintf('Rat%s', rat_id);
    if ~isfield(SR_DATA, ratFld); continue; end
    rat_entry = SR_DATA.(ratFld);

    for sIdx = 1:numel(sessions)
        sess_name = sessions{sIdx};
        for rg = 1:numel(regions)
            region = regions{rg};
            for wn = 1:numel(windows)
                win_fld = windows{wn};
                row = row + 1;

                Rat(row)     = rat_id;
                Age(row)     = age;
                Session(row) = sess_name;
                Region(row)  = region;
                Window(row)  = win_lbls{wn};

                % L2 and L5
                [out_L2, det_win_L2] = compute_for_layer(rat_entry, ...
                    rat_entry.L2_shank, sess_name, region, win_fld, T_VEC);
                [out_L5, det_win_L5] = compute_for_layer(rat_entry, ...
                    rat_entry.L5_shank, sess_name, region, win_fld, T_VEC);

                if out_L2.valid
                    Peak_Latency_L2_ms(row)   = out_L2.peak_lat_ms;
                    Peak_Amplitude_L2_mV(row) = out_L2.peak_amp_mV;
                end
                if out_L5.valid
                    Peak_Latency_L5_ms(row)   = out_L5.peak_lat_ms;
                    Peak_Amplitude_L5_mV(row) = out_L5.peak_amp_mV;
                end

                % Union across layers (skip layer if invalid)
                [hm_s, hm_e] = union_crossings(out_L2, out_L5, ...
                    'halfmax_start_ms', 'halfmax_end_ms');
                [zr_s, zr_e] = union_crossings(out_L2, out_L5, ...
                    'zero_start_ms',    'zero_end_ms');

                HalfMax_Start_ms(row)    = hm_s;
                HalfMax_End_ms(row)      = hm_e;
                HalfMax_Duration_ms(row) = hm_e - hm_s;
                Zero_Start_ms(row)       = zr_s;
                Zero_End_ms(row)         = zr_e;
                Zero_Duration_ms(row)    = zr_e - zr_s;
            end
        end
    end
end

% Trim unused rows
Rat(row+1:end) = []; Age(row+1:end) = []; Session(row+1:end) = [];
Region(row+1:end) = []; Window(row+1:end) = [];
Peak_Latency_L2_ms(row+1:end) = []; Peak_Amplitude_L2_mV(row+1:end) = [];
Peak_Latency_L5_ms(row+1:end) = []; Peak_Amplitude_L5_mV(row+1:end) = [];
HalfMax_Start_ms(row+1:end) = []; HalfMax_End_ms(row+1:end) = [];
HalfMax_Duration_ms(row+1:end) = [];
Zero_Start_ms(row+1:end) = []; Zero_End_ms(row+1:end) = [];
Zero_Duration_ms(row+1:end) = [];

T = table(Rat, Age, Session, Region, Window, ...
    Peak_Latency_L2_ms, Peak_Amplitude_L2_mV, ...
    Peak_Latency_L5_ms, Peak_Amplitude_L5_mV, ...
    HalfMax_Start_ms, HalfMax_End_ms, HalfMax_Duration_ms, ...
    Zero_Start_ms, Zero_End_ms, Zero_Duration_ms);

writetable(T, save_path);
fprintf('Peak window times CSV: %s  (%d rows)\n', save_path, height(T));

end


% =========================================================================
f
