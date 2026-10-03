
% anglesDeg = [0 45 90 135 180 225 270 315];
% Lp = [matrix: rows = angels, columns = frequencies]

function [D_source, L_ref] = calc_Dsource(anglesDeg, Lp_angles)

    anglesDeg = anglesDeg(:);

    if size(Lp_angles,1) ~= length(anglesDeg)
        error('Number of rows in Lp must match number of angles.');
    end

    linearLevels = 10.^(Lp_angles ./ 10);

    L_ref = 10 .* log10(mean(linearLevels,1));

    D_source = Lp_angles - L_ref;

end