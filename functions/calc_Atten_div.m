% Geometry attenuation for a single-point sound source in a free space, 
% assuming a spherical propogation. For a hemi-spherical, use the relevant
% varialbe in the directivity attenuation function.
function A_div = calc_Atten_div(d)
    
    d0 = 1; % [m], constant 
    S0 = 4*pi*d0^2; % surface area of a sphere

    % Lp = 10*log10((W/I0)/(4*pi*d^2))
    %       = 10*log10(W/I0) - 10*log10(4*pi*d^2)
    %       = L_W - 10*log10(4*pi) - 20*log10(d)
    %       = L_W - 11 - 20*log10(d)

    % Therefore the attenuation is:
    A_div = 20*log10(d/d0) + 11;

    % bigger d -> bigger A_div
end