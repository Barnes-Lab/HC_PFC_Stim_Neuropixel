function ds = downsampleLFP(dataVolts, Fs, Fs_ds)
%  Downsample to Fs_ds then zero-phase band-pass (1–300 Hz)
% applying padded filtering to minimize edge artifacts.
% Inputs:
%   dataVolts: (nChan x nTime) raw voltage data per trial(Volts)
%   Fs:original sampling rate (Hz)
%   Fs_ds: downsample rate (Hz)
% Output:
%   dataFiltered  (nChan x nTime) 

   
    % Downsample factors
    [P, Q] = rat(Fs_ds / Fs, 1e-8);

    % 2) Downsample with anti-aliasing filter function
    ds = resample(dataVolts.', P, Q);  
    ds = ds.';  
    % nT2 = size(ds, 2);
    % 
    % % 3) Determine padding length: use smaller of Fs_ds or nT2-1
    % padSamples = min(Fs_ds, nT2-1);
    % if padSamples < 1
    %     padSamples = 0;  % no padding if too short
    % end
    % 
    % % 4) Pad on both ends by mirroring to reduce transients
    % if padSamples > 0
    %     padLeft  = ds(:, padSamples:-1:1);
    %     padRight = ds(:, end:-1:end-padSamples+1);
    %     dsPad = [padLeft, ds, padRight];  % [nCh x (nT2 + 2*padSamples)]
    % else
    %     dsPad = ds;
    % end
    % 
    % % 5) Design 12th-order Butterworth band-pass 0.5–300 Hz
    % Wn = [0.5 300] / (Fs_ds / 2);
    % [b, a] = butter(12, 300/(Fs_ds / 2));
    % 
    % % 6) Zero-phase filter padded data
    % dsFiltPad = filtfilt(b, a, dsPad.').';  % [nCh x padded length]
    % 
    % % 7) Remove padding to recover original downsampled length
    % if padSamples > 0
    %     dataFilt = dsFiltPad(:, padSamples+1 : padSamples + nT2);
    % else
    %     dataFilt = dsFiltPad;
    % end
end
