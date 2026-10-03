function sigBand = octaveFilter_signalToolbox(sig,fs,f_center)

    G = 10^(3/10);

    f_low  = f_center * G^(-1/2);
    f_high = f_center * G^(1/2);

    if f_high >= fs/2
        error('Sampling rate is too low for this octave band.');
    end

    sig = sig(:);

    Wn = [f_low f_high]/(fs/2);

    filterOrder = 4;

    [z,p,k] = butter(filterOrder,Wn,'bandpass');

    sos = zp2sos(z,p,k);

    sigBand = sosfilt(sos,sig);

end