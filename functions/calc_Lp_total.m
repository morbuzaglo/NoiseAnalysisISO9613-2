function Lp_total = calc_Lp_total(sig,fs,f,t1,t2,bandType)

    Lp_total = zeros(size(f));

    for i = 1:length(f)
        Lp_total(i) = calc_Lp_band(sig,fs,f(i),t1,t2,bandType);
    end

end