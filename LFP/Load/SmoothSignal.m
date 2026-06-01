function dataSmooth = SmoothSignal(dataFilt, Fs)
    % Smooth data either through a butterworth filter or through sgolay filter
    
    % % Get no. of samples
    nT2 = size(dataFilt, 2);

    %  Fs or nT2-1
    padSamples = min(Fs, nT2-1);
    if padSamples < 1
        padSamples = 0;  % no padding if too short
    end

    % Pad on both ends 
    if padSamples > 0
        padLeft  = dataFilt(:, padSamples:-1:1);
        padRight = dataFilt(:, end:-1:end-padSamples+1);
        dsPad = [padLeft, dataFilt, padRight];  
    else
        dsPad = dataFilt;
    end

    % 12th-order Butterworth low pass filter 300 Hz
    Fc       = 300; order = 12;
    [b, a]   = butter(order, Fc/(Fs/2), 'low');

    %  Zero-phase filter padded data
    dsFiltPad = filtfilt(b, a, dsPad.').'; 

    % Remove padding 
    if padSamples > 0
        dataSmooth = dsFiltPad(:, padSamples+1 : padSamples + nT2);
    else
        dataSmooth = dsFiltPad;
    end

    % [nCh, ~] = size(dataFilt);
    % dataSmooth=zeros(size(dataFilt));
    % % 1)Butter worth filter 
    % Fc       = 300; order = 12;
    % [b, a]   = butter(order, Fc/(Fs/2), 'low');
    % dataSmooth = zeros(size(dataFilt));
    % for ch = 1:nCh
    %     dataSmooth(ch, :) = filtfilt(b, a, dataFilt(ch, :));
    % 
    % end

    % 2) Savitzky Golay
    % polyOrder   = 3;
    % frameLength = round(0.003 * Fs)+1;  % 3 ms window-Odd
    % if mod(frameLength,2)==0
    %     frameLength = frameLength+1;
    % end
    % for ch = 1:nCh
    %     dataSmooth(ch, :) =sgolayfilt(dataFilt(ch, :), polyOrder, frameLength);
    % end
end
