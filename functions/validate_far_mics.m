% Validate and check all of the far microphones at the experiment, in order
% to see if the prediction is good or if adjuments/further calibration are
% required in the predeiction model. 

function results = validate_far_mics( ...
    f,L_W,calibAnglesDeg,D_source, ...
    Lp_total,Lp_bg,thetaDeg, ...
    D_surface,A_div,A_atm,A_gr,A_bar)

    nMics = size(Lp_total,1);
    nFreq = length(f);

    if size(Lp_total,2) ~= nFreq
        error('Lp_total must have one column per frequency band.');
    end

    if ~isequal(size(Lp_bg),size(Lp_total))
        error('Lp_bg must have the same size as Lp_total.');
    end

    thetaDeg = thetaDeg(:);

    if length(thetaDeg) ~= nMics
        error('thetaDeg must contain one angle per microphone.');
    end

    if length(L_W) ~= nFreq
        error('L_W must contain one value per frequency band.');
    end

    if isscalar(D_surface)
        D_surface = D_surface .* ones(nMics,nFreq);
    elseif isvector(D_surface) && length(D_surface) == nFreq
        D_surface = repmat(D_surface(:).',nMics,1);
    elseif ~isequal(size(D_surface),[nMics nFreq])
        error('D_surface must be scalar, 1 x nFreq, or nMics x nFreq.');
    end

    if isscalar(A_div)
        A_div = A_div .* ones(nMics,nFreq);
    elseif isvector(A_div) && length(A_div) == nMics
        A_div = repmat(A_div(:),1,nFreq);
    elseif ~isequal(size(A_div),[nMics nFreq])
        error('A_div must be scalar, nMics x 1, or nMics x nFreq.');
    end

    if isscalar(A_atm)
        A_atm = A_atm .* ones(nMics,nFreq);
    elseif isvector(A_atm) && length(A_atm) == nFreq
        A_atm = repmat(A_atm(:).',nMics,1);
    elseif ~isequal(size(A_atm),[nMics nFreq])
        error('A_atm must be scalar, 1 x nFreq, or nMics x nFreq.');
    end

    if isscalar(A_gr)
        A_gr = A_gr .* ones(nMics,nFreq);
    elseif isvector(A_gr) && length(A_gr) == nFreq
        A_gr = repmat(A_gr(:).',nMics,1);
    elseif ~isequal(size(A_gr),[nMics nFreq])
        error('A_gr must be scalar, 1 x nFreq, or nMics x nFreq.');
    end

    if isscalar(A_bar)
        A_bar = A_bar .* ones(nMics,nFreq);
    elseif isvector(A_bar) && length(A_bar) == nFreq
        A_bar = repmat(A_bar(:).',nMics,1);
    elseif ~isequal(size(A_bar),[nMics nFreq])
        error('A_bar must be scalar, 1 x nFreq, or nMics x nFreq.');
    end

    Lp_measured = nan(nMics,nFreq);
    validBG = false(nMics,nFreq);

    Lp_predicted = nan(nMics,nFreq);
    D_source_est = nan(nMics,nFreq);

    for i = 1:nMics

        [Lp_measured(i,:),validBG(i,:)] = ...
            remove_background( ...
                Lp_total(i,:), ...
                Lp_bg(i,:));

        D_source_est(i,:) = ...
            estimate_Dsource_at_angle( ...
                calibAnglesDeg, ...
                D_source, ...
                thetaDeg(i));

        Lp_predicted(i,:) = ...
            L_W + ...
            D_source_est(i,:) + ...
            D_surface(i,:) - ...
            A_div(i,:) - ...
            A_atm(i,:) - ...
            A_gr(i,:) - ...
            A_bar(i,:);

    end

    Delta_L = Lp_measured - Lp_predicted;

    valid = ...
        validBG & ...
        ~isnan(Lp_measured) & ...
        ~isnan(Lp_predicted);

    Bias = nan(nMics,1);
    MAE = nan(nMics,1);
    RMSE = nan(nMics,1);

    for i = 1:nMics

        idx = valid(i,:);

        if any(idx)

            Bias(i) = mean( ...
                Delta_L(i,idx));

            MAE(i) = mean( ...
                abs(Delta_L(i,idx)));

            RMSE(i) = sqrt( ...
                mean(Delta_L(i,idx).^2));

        end

    end

    Bias_freq = nan(1,nFreq);
    MAE_freq = nan(1,nFreq);
    RMSE_freq = nan(1,nFreq);

    for j = 1:nFreq

        idx = valid(:,j);

        if any(idx)

            Bias_freq(j) = mean( ...
                Delta_L(idx,j));

            MAE_freq(j) = mean( ...
                abs(Delta_L(idx,j)));

            RMSE_freq(j) = sqrt( ...
                mean(Delta_L(idx,j).^2));

        end

    end

    A_weight = ...
        [-26.2 -16.1 -8.6 -3.2 0 1.2 1.0 -1.1];

    if nFreq ~= length(A_weight)
        error('A-weighting vector does not match frequency vector.');
    end

    LpA_measured = nan(nMics,1);
    LpA_predicted = nan(nMics,1);
    Delta_LA = nan(nMics,1);

    for i = 1:nMics

        idx = valid(i,:);

        if any(idx)

            LpA_measured(i) = ...
                10*log10(sum( ...
                    10.^(( ...
                    Lp_measured(i,idx) + ...
                    A_weight(idx))/10)));

            LpA_predicted(i) = ...
                10*log10(sum( ...
                    10.^(( ...
                    Lp_predicted(i,idx) + ...
                    A_weight(idx))/10)));

            Delta_LA(i) = ...
                LpA_measured(i) - ...
                LpA_predicted(i);

        end

    end

    results.f = f;

    results.Lp_measured = Lp_measured;
    results.Lp_predicted = Lp_predicted;
    results.Delta_L = Delta_L;

    results.validBG = validBG;
    results.valid = valid;

    results.D_source_est = D_source_est;

    results.D_surface = D_surface;
    results.A_div = A_div;
    results.A_atm = A_atm;
    results.A_gr = A_gr;
    results.A_bar = A_bar;

    results.Bias = Bias;
    results.MAE = MAE;
    results.RMSE = RMSE;

    results.Bias_freq = Bias_freq;
    results.MAE_freq = MAE_freq;
    results.RMSE_freq = RMSE_freq;

    results.LpA_measured = LpA_measured;
    results.LpA_predicted = LpA_predicted;
    results.Delta_LA = Delta_LA;

end