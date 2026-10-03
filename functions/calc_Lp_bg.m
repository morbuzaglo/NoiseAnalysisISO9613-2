function Lp_bg = calc_Lp_bg(sig_bg,fs,f,t1,t2,bandType)

    Lp_bg = zeros(size(f));

    for i = 1:length(f)
        Lp_bg(i) = calc_Lp_band(sig_bg,fs,f(i),t1,t2,bandType);
    end

end