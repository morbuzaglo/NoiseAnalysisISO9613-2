function sigBand = octaveFilter_audioToolbox(sig,fs,f_center)

    if exist('octaveFilter','class') ~= 8 && ...
       exist('octaveFilter','file') ~= 2
        error('Audio Toolbox octaveFilter is not available.');
    end

    G = 10^(3/10);

    f_high = f_center * G^(1/2);

    if f_high >= fs/2
        error('Sampling rate is too low for this octave band.');
    end

    octFilt = octaveFilter( ...
        'CenterFrequency',f_center, ...
        'Bandwidth','1 octave', ...
        'SampleRate',fs, ...
        'FilterOrder',8, ...
        'Oversample',true);

    sigBand = octFilt(sig(:));

end