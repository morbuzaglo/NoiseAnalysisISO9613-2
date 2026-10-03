% THIS VALUE ADDS TO THE ATTENUATION, and CAN be used in the case of not
% downwind or the wind soeed in the other way sometimes.. something like
% that. don't believe me? GO READ IT.

% NOT as a part of the (frequency-dependent) attenuations addition, but use
% it ONLY at the end AFTER comuting the equivalent final A-weighted sound
% level value.

% C0 - around the values: 0 dB -> 5 dB, most of the time. must be computed
% using a meteorological calculation.

function C_met = calc_Cmet(h_s, h_r, d_p, C0)

    % Compute the long meteorological correction value:
    if d_p <= 10*(h_s + h_r) % Too short of a distance for the meteorological conditions to affect
        C_met = 0;
    else
        C_met = C0*(1 - 10*(h_s + h_r)/d_p);
    end

end
