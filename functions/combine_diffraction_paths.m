% For use in the barrier attenuation function.

% Inputs format:
%   - A_paths = [
%               A_bar_top
%               A_bar_left
%               A_bar_right
%                         ]; column = frequency, row = paths

function A_bar_total = combine_diffraction_paths(A_paths)

    validPaths = ~isnan(A_paths);

    transmission = zeros(size(A_paths));

    transmission(validPaths) = ...
        10.^(-A_paths(validPaths) ./ 10);

    totalTransmission = sum(transmission, 1);

    A_bar_total = -10 .* log10(totalTransmission);

    A_bar_total = max(A_bar_total, 0);

end