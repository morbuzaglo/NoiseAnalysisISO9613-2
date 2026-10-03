function [Lp_source,valid] = remove_background(Lp_total,Lp_bg)

    delta = Lp_total - Lp_bg;

    valid = delta >= 3;

    Lp_source = nan(size(Lp_total));

    sourceEnergy = 10.^(Lp_total(valid)/10) - 10.^(Lp_bg(valid)/10);

    Lp_source(valid) = 10*log10(sourceEnergy);

end