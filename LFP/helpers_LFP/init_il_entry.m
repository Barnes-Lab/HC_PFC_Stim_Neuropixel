function il = init_il_entry(n_stim, n_chan, mono_window_ms)
% Initialise an IL session entry with empty mono peak struct
% SS 2026

il = struct();
il.mono  = empty_il_mono(n_stim, n_chan, mono_window_ms);
il.strong_stims  = [];
il.excluded_channels = [];
il.channel_mode  = 'none';
il.frac_walled  = NaN;
end
