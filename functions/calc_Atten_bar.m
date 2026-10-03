% Thumb rules:
%   - Higher frequincies -> bigger A_bar
%   - Lower frequencies -> easily bends around obstacles
%   - More screening (less Line-Of-Sight) -> bigger A_bar

% Important:
%   - A_bar INCLUDES A_gr, unless A_gr <=0 (applifies signal) or in the 
%     case of vertical edges. In those cases, both should be added to
%     the total attenuation: A_diff + A_gr.

% Inputs units:
%   - frequency f [Hz]
%   - speed of sound c [m/s] @20 degrees C
%   - difference between diffracted route to straight route [m]
%   - The distance between the first and last barriers' edges [m]
%   - diffractionType = 'top' or 'lateral'

function [A_bar, D_z, K_met, C3, z_min] = calc_Atten_bar(f, c, z, e, d_ss, d_sr, d, A_gr, diffractionType)

    lambda = c ./ f; % Wave length
    C2 = 20; % constant (sometimes = 40, not this case)

    if e == 0 % if only one edge exist
        C3 = ones(size(f));
    else % more than one barrier edge
        C3 = (1 + (5 .* lambda ./ e).^2) ./ ...
             (1/3 + (5 .* lambda ./ e).^2);
    end

    z_min = -2 .* lambda ./ (C2 .* C3);
    K_met = ones(size(f)); % meteorological correction factor

    valid = z > z_min;

    if strcmpi(diffractionType, 'top') % above the barrier
        K_met(valid) = exp( ...
            -(1/2000) .* sqrt( ...
            ((max(d_ss, d_sr) + e) .* min(d_ss, d_sr) .* d) ./ ...
            (2 .* (z - z_min(valid))) ...
            ) ...
            );
    elseif strcmpi(diffractionType, 'lateral') % right/left of the barrier
        K_met(:) = 1;
    else
        error('diffractionType must be ''top'' or ''lateral''.');
    end

    % Diffraction attenuation:
    D_z = zeros(size(f));
    D_z(valid) = 10 .* log10( ...
        1 + ...
        (2 + C2 ./ lambda(valid)) .* ...
        C3(valid) .* ...
        z .* ...
        K_met(valid) ...
        );

    if isscalar(A_gr)
        A_gr = A_gr .* ones(size(f));
    end

    if strcmpi(diffractionType, 'top')
        A_bar = zeros(size(f));

        positiveGround = A_gr > 0;

        A_bar(positiveGround) = D_z(positiveGround) - A_gr(positiveGround);
        A_bar(~positiveGround) = D_z(~positiveGround);

        A_bar = max(A_bar, 0);
    else
        A_bar = max(D_z, 0); % is D_z<0 -> we must take 0 (by ISO)
    end

end