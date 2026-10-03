function [Lp,p_rms,sigBand,methodUsed] = calc_Lp_band(sig,fs,f_center,t1,t2,bandType)

    p0 = 20e-6;

    i1 = max(1,round(t1*fs)+1);
    i2 = min(length(sig),round(t2*fs));

    if i2 <= i1
        error('Invalid time interval.');
    end

    sigSegment = sig(i1:i2);

    switch lower(bandType)

        case 'octave'

            G = 10^(3/10);

            f_low  = f_center*G^(-1/2);
            f_high = f_center*G^(1/2);

        otherwise
            error('This implementation currently supports octave bands only.');

    end

    if f_high >= fs/2
        error('Sampling rate is too low for this frequency band.');
    end

    if fs < 2.5*f_high
        warning('Sampling rate has limited margin above the upper band frequency.');
    end

    T = length(sigSegment)/fs;

    if T*f_low < 20
        warning('Time window contains fewer than 20 cycles of the lower band edge.');
    end

    % Filter using a function I wrote:
    [sigBand,methodUsed] = octaveFilter_auto(sigSegment,fs,f_center);
    p_rms = sqrt(mean(sigBand.^2));

    Lp = 20*log10(p_rms/p0);

end