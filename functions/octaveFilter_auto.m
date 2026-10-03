% Choses the method to use based on what toolboxes I have on my computer.

function [sigBand,methodUsed] = octaveFilter_auto(sig,fs,f_center)

    hasAudioToolbox = ...
        exist('octaveFilter','class') == 8 || ...
        exist('octaveFilter','file') == 2;

    hasSignalToolbox = ...
        exist('butter','file') == 2 && ...
        exist('sosfiltfilt','file') == 2;

    if hasAudioToolbox

        sigBand = octaveFilter_audioToolbox(sig,fs,f_center);
        methodUsed = 'Audio Toolbox octaveFilter';

    elseif hasSignalToolbox

        sigBand = octaveFilter_signalToolbox_zeroPhase(sig,fs,f_center);
        methodUsed = 'Signal Processing Toolbox Butterworth';

    else

        sigBand = octaveFilter_IEC(sig,fs,f_center);
        methodUsed = 'Custom FIR implementation';

    end

end