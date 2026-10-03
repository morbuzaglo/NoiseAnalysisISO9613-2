% Predict sound level at a distant microphone, using the data extracted
% from the far ones.
% Attenuation values are computed before.
function [Lp_pred,D_source_est] = predict_Lp_at_mic( ...
    f,L_W,thetaDeg,calibAnglesDeg,D_source, ...
    D_surface,A_div,A_atm,A_gr,A_bar)

    nFreq = length(f);

    if isscalar(D_surface)
        D_surface = D_surface .* ones(1,nFreq);
    end

    if isscalar(A_div)
        A_div = A_div .* ones(1,nFreq);
    end

    if isscalar(A_atm)
        A_atm = A_atm .* ones(1,nFreq);
    end

    if isscalar(A_gr)
        A_gr = A_gr .* ones(1,nFreq);
    end

    if isscalar(A_bar)
        A_bar = A_bar .* ones(1,nFreq);
    end

    D_source_est = estimate_Dsource_at_angle( ...
        calibAnglesDeg,D_source,thetaDeg);

    Lp_pred = L_W + D_source_est + D_surface - A_div - A_atm - A_gr - A_bar;

end
