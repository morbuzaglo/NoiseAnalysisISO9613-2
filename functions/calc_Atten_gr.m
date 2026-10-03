% s - source
% r - reciever
% 0 < G < 1, where [0 = soft ground] and [1 = Hard/reflecting ground]

% ONLY relevant for frequencies: 
% [63, 125, 250, 500, 1000, 2000, 4000, 8000] Hz

% Thumb rules:
%   - if the source and reciever are very far from each other,
%     meaning dp >> hs, hr -> bigger A_gr
  
function [A_gr, A_s, A_r, A_m] = calc_Atten_gr(f, h_s, h_r, d_p, G_s, G_m, G_r)
    % q = a factor that decides how much the middle area between the source
    % and the reciever really contributes:
    if d_p <= 30*(h_s + h_r)
        q = 0; % No middle area
    else
        q = 1 - 30*(h_s + h_r)/d_p; % Some middle area exist
    end
    
    % Help expressions - source side:
    a_s = 1.5 + ...
          3.0*exp(-0.12*(h_s - 5)^2)*(1 - exp(-d_p/50)) + ...
          5.7*exp(-0.09*h_s^2)*(1 - exp(-2.8e-6*d_p^2));
    b_s = 1.5 + ...
          8.6*exp(-0.09*h_s^2)*(1 - exp(-d_p/50));
    c_s = 1.5 + ...
          14.0*exp(-0.46*h_s^2)*(1 - exp(-d_p/50));
    d_s = 1.5 + ...
          5.0*exp(-0.9*h_s^2)*(1 - exp(-d_p/50));

    % Help expressions - reciever side:
    a_r = 1.5 + ...
          3.0*exp(-0.12*(h_r - 5)^2)*(1 - exp(-d_p/50)) + ...
          5.7*exp(-0.09*h_r^2)*(1 - exp(-2.8e-6*d_p^2));
    b_r = 1.5 + ...
          8.6*exp(-0.09*h_r^2)*(1 - exp(-d_p/50));
    c_r = 1.5 + ...
          14.0*exp(-0.46*h_r^2)*(1 - exp(-d_p/50));
    d_r = 1.5 + ...
          5.0*exp(-0.9*h_r^2)*(1 - exp(-d_p/50));
    
    A_s = zeros(size(f)); % Source ground side attenuation contribution
    A_r = zeros(size(f)); % Reciever ground side attenuation contribution
    A_m = zeros(size(f)); % Middle ground side attenuation contribution
    
    for i = 1:numel(f)
    
        switch f(i)
            case 63
                A_s(i) = -1.5;
                A_r(i) = -1.5;
                A_m(i) = -3*q;
            case 125
                A_s(i) = -1.5 + G_s*a_s;
                A_r(i) = -1.5 + G_r*a_r;
                A_m(i) = -3*q*(1 - G_m);
            case 250
                A_s(i) = -1.5 + G_s*b_s;
                A_r(i) = -1.5 + G_r*b_r;
                A_m(i) = -3*q*(1 - G_m);
            case 500
                A_s(i) = -1.5 + G_s*c_s;
                A_r(i) = -1.5 + G_r*c_r;
                A_m(i) = -3*q*(1 - G_m);
            case 1000
                A_s(i) = -1.5 + G_s*d_s;
                A_r(i) = -1.5 + G_r*d_r;
                A_m(i) = -3*q*(1 - G_m);
            case {2000, 4000, 8000}
                A_s(i) = -1.5*(1 - G_s);
                A_r(i) = -1.5*(1 - G_r);
                A_m(i) = -3*q*(1 - G_m);
            otherwise
                error('Frequency must be one of: 63, 125, 250, 500, 1000, 2000, 4000, 8000 Hz.')
        end
    end

    % Total ground attenuation, NOT final (correction K_geo needed):
    A_gr_prime = A_s + A_r + A_m; 
    
    % Correction factor for geometry:
    K_geo = (d_p^2 + (h_s - h_r)^2) / ...
            (d_p^2 + (h_s + h_r)^2);
    
    % Total (and final) ground attenuation:
    A_gr = -10*log10(1 + (10.^(-A_gr_prime/10) - 1)*K_geo);

end