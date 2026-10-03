% 'normal' incidence is more restrictive then 'diffuse', so better use it.

% Thumb rules:
%   - Higher density or thickness -> higher mass attenuation R
%   - Higher frequency -> higher mass attenuation R

function R = calc_Atten_masslaw(f,rho,t,incidence)

    m_prime = rho * t; % [kg/m^2]

    % C - a factor related to the incidence angles with the mass:
    switch lower(incidence)
        case 'normal'
            C = 42;
        case 'diffuse'
            C = 47;
        otherwise
            error('incidence must be ''normal'' or ''diffuse''.');
    end

    % Mass-law attenuation (frequency-dependent) value:
    R = 20 .* log10(m_prime .* f) - C;

end