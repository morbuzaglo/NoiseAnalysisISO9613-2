% This atmospheric attenuation DOES NOT indlude:
%   - Wind refraction
%   - Temperature-gradient refraction
%   - Turbulence
%   - Ground effects

% Relevant for:
%   - Tempratures: -20C < T < 50C
%   - Relative humidity: 0.1 < RH < 1
%   - Sound frequencies: 50Hz < f < 10,000Hz

% Inputs units:
%   - f [Hz]
%   - T [degrees C]
%   - RH [%] -> 0<RH<1
%   - p [Pa]
%   - d [m]

% Thumb rules:
%   - Higher sound frequency -> bigger A_atm (STRONG)
%   - bigger distance -> bigger A_atm (linear)
%   - higher temprature -> inconclusive (mostly bigger A_atm)
%   - higher RH -> inconclusive

function [A_atm, alpha] = calc_Atten_atm(f, T, RH, p, d)

    T_K = T + 273.15; % [K]
    
    T0 = 293.15; % Reference temp value [K]
    p0 = 101325; % Reference pressure value [Pa]
    
    h = RH .* 10.^(-6.8346 .* (273.16 ./ T_K).^1.261 + 4.6151) .* (p0 ./ p);
    
    % Relaxation frequency of Oxigen:
    frO = (p ./ p0) .* ...
          (24 + 4.04e4 .* h .* (0.02 + h) ./ (0.391 + h));
    
    % Relaxation frequency of Nitrogen:
    frN = (p ./ p0) .* (T_K ./ T0).^(-0.5) .* ...
          (9 + 280 .* h .* exp(-4.17 .* ((T_K ./ T0).^(-1/3) - 1)));
    
    % The atmospheric attenuation factor (frequency-dependent):
    alpha = 1000 .* 8.686 .* f.^2 .* ( ...
            1.84e-11 .* (p0 ./ p) .* sqrt(T_K ./ T0) + ...
            (T_K ./ T0).^(-2.5) .* ...
            ( ...
            0.01275 .* exp(-2239.1 ./ T_K) ./ (frO + f.^2 ./ frO) + ...
            0.1068 .* exp(-3352 ./ T_K) ./ (frN + f.^2 ./ frN) ...
            ) ...
            );
    
    % Atmospheric attenuation:
    A_atm = alpha .* d ./ 1000;

end