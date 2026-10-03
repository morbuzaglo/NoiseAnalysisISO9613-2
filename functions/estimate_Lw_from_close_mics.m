% Estimate the LW power level [dB] and source directivity from the close
% mics in the experiment, assuming no barrier between them.

function [L_W,L_W_std,L_W_mics,D_source,Lp_source,validBG] = ...
    estimate_Lw_from_close_mics( ...
    f,Lp_total,Lp_bg,d,d_p,h_s,h_r,G_s,G_m,G_r,T_C,RH,p,D_surface,anglesDeg)

    nMics = size(Lp_total,1);
    nFreq = length(f);

    anglesDeg = anglesDeg(:);

    if size(Lp_total,2) ~= nFreq || size(Lp_bg,2) ~= nFreq
        error('Lp_total and Lp_bg must have one column per frequency band.');
    end

    if size(Lp_bg,1) ~= nMics
        error('Lp_total and Lp_bg must contain the same number of microphones.');
    end

    if length(anglesDeg) ~= nMics
        error('anglesDeg must contain one angle for each microphone.');
    end

    if isscalar(D_surface)
        D_surface = D_surface .* ones(nMics,nFreq);
    elseif isvector(D_surface) && length(D_surface) == nFreq
        D_surface = repmat(D_surface(:).',nMics,1);
    end

    if isscalar(h_r)
        h_r = h_r .* ones(nMics,1);
    end

    if isscalar(G_s)
        G_s = G_s .* ones(nMics,1);
    end

    if isscalar(G_m)
        G_m = G_m .* ones(nMics,1);
    end

    if isscalar(G_r)
        G_r = G_r .* ones(nMics,1);
    end

    deltaBG = Lp_total - Lp_bg;
    validBG = deltaBG >= 3;

    Lp_source = nan(size(Lp_total));

    for i = 1:nMics
        for j = 1:nFreq

            if validBG(i,j)

                sourceEnergy = 10.^(Lp_total(i,j)/10) - ...
                               10.^(Lp_bg(i,j)/10);

                if sourceEnergy > 0
                    Lp_source(i,j) = 10*log10(sourceEnergy);
                end

            end
        end
    end

    L_W_mics = nan(nMics,nFreq);

    for i = 1:nMics

        A_div = calc_Atten_div(d(i));

        [A_atm,~] = calc_Atten_atm(f,T_C,RH,p,d(i));

        [A_gr,~,~,~] = calc_Atten_gr( ...
            f,h_s,h_r(i),d_p(i),G_s(i),G_m(i),G_r(i));

        L_W_mics(i,:) = Lp_source(i,:) - D_surface(i,:) + ...
                        A_div + A_atm + A_gr;

    end

    L_W = nan(1,nFreq);

    for j = 1:nFreq

        valid = ~isnan(L_W_mics(:,j));

        if any(valid)
            L_W(j) = 10*log10(mean(10.^(L_W_mics(valid,j)/10)));
        end

    end

    D_source = L_W_mics - L_W;

    L_W_std = std(L_W_mics,0,1,'omitnan');

    [anglesDeg,order] = sort(anglesDeg);
    D_source = D_source(order,:);
    L_W_mics = L_W_mics(order,:);
    Lp_source = Lp_source(order,:);
    validBG = validBG(order,:);

end