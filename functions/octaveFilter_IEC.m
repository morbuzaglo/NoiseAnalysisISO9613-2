% A function we wrote in order to avoid using the toolboxes.
% according to IEC 61260-1, class 1.

function sigBand = octaveFilter_IEC(sig,fs,f_center)

    G = 10^(3/10);

    f_low  = f_center * G^(-1/2);
    f_high = f_center * G^(1/2);

    if f_high >= fs/2
        error('Sampling rate is too low for this octave band.');
    end

    sig = sig(:);

    transitionFraction = 0.12;

    df_low  = transitionFraction*f_low;
    df_high = transitionFraction*f_high;

    df = min(df_low,df_high);

    Astop = 80;

    Dw = 2*pi*df/fs;

    N = ceil((Astop-8)/(2.285*Dw));

    if mod(N,2) ~= 0
        N = N + 1;
    end

    n = 0:N;
    M = N/2;

    fc1 = f_low/fs;
    fc2 = f_high/fs;

    h = 2*fc2*sinc_custom(2*fc2*(n-M)) - ...
        2*fc1*sinc_custom(2*fc1*(n-M));

    beta = kaiser_beta(Astop);

    x = 2*n/N - 1;

    w = besseli(0,beta*sqrt(1-x.^2)) / besseli(0,beta);

    h = h .* w;

    Hc = sum(h .* exp(-1i*2*pi*f_center/fs*n));

    h = h / abs(Hc);

    sigBand = conv(sig,h(:),'same');

end


function y = sinc_custom(x)

    y = ones(size(x));

    idx = abs(x) > eps;

    y(idx) = sin(pi*x(idx))./(pi*x(idx));

end


function beta = kaiser_beta(A)

    if A > 50
        beta = 0.1102*(A-8.7);

    elseif A >= 21
        beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);

    else
        beta = 0;
    end

end