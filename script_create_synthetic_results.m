clear
clc
close all

%% ============================================================
%  SYNTHETIC EXPERIMENT GENERATOR
% =============================================================

outputFolder = fullfile(pwd,'experiment');

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%% General parameters

fs = 48000;
duration = 20;

t = (0:1/fs:duration-1/fs).';

f = [63 125 250 500 1000 2000 4000 8000];

T_C = 20;
RH = 70;
p = 101325;

h_s = 1.5;

%% Microphone geometry

distances = ...
    [7 14 21 50 250 7 14 21 50 7 14 21 50];

heights = ...
    [1.5 1.5 1.5 1.5 -11 ...
     1.5 1.5 1.5 1.5 ...
     1.5 1.5 1.5 1.5];

angles = ...
    [0 0 0 0 0 ...
     120 120 120 120 ...
     240 240 240 240];

nMics = length(distances);

%% ============================================================
%  TRUE SOURCE CHARACTERISTICS
% =============================================================

L_W_true = ...
    [94 98 101 103 104 102 98 92];

%% Source directivity

D0 = ...
    [ 2  2  3  4  5  5  4  3];

D120 = ...
    [-1  0  0 -1 -2 -2 -3 -3];

D240 = ...
    [-2 -2 -3 -3 -3 -2 -1  0];

%% ============================================================
%  SIGNAL ACTIVITY
% =============================================================

sourceStart = 5;
sourceEnd = 15;

sourceMask = ...
    t >= sourceStart & t <= sourceEnd;

%% ============================================================
%  BACKGROUND SPECTRUM
% =============================================================

% Typical quiet residential-ish background

Lp_bg = ...
    [46 41 36 31 27 23 19 16];

%% ============================================================
%  GENERATE EACH MICROPHONE
% =============================================================

rng(1)

for i = 1:nMics

    d_p = distances(i);

    h_r = heights(i);

    d = sqrt( ...
        d_p^2 + ...
        (h_s-h_r)^2);

    %% Directivity

    if angles(i) == 0
        D_source = D0;

    elseif angles(i) == 120
        D_source = D120;

    else
        D_source = D240;
    end

    %% Propagation attenuation

    A_div = calc_Atten_div(d);

    [A_atm,~] = ...
        calc_Atten_atm( ...
            f,T_C,RH,p,d);

    [A_gr,~,~,~] = ...
        calc_Atten_gr( ...
            f, ...
            h_s, ...
            h_r, ...
            d_p, ...
            1,1,1);

    %% Extra topography effect for Mic 5

    A_bar = zeros(size(f));

    if i == 5

        A_bar = ...
            [0 0 1 3 6 9 12 15];

    end

    %% SPL at microphone

    Lp_source = ...
        L_W_true + ...
        D_source - ...
        A_div - ...
        A_atm - ...
        A_gr - ...
        A_bar;

    %% ========================================================
    %  BUILD TIME-DOMAIN SIGNAL
    % =========================================================

    sigSource = zeros(size(t));

    sigBackground = zeros(size(t));

    for j = 1:length(f)

        %% Source amplitude

        p_rms_source = ...
            20e-6 * ...
            10^(Lp_source(j)/20);

        p_peak_source = ...
            sqrt(2)*p_rms_source;

        phaseSource = ...
            2*pi*rand;

        toneSource = ...
            p_peak_source * ...
            sin( ...
                2*pi*f(j)*t + ...
                phaseSource);

        sigSource = ...
            sigSource + ...
            toneSource .* sourceMask;

        %% Background amplitude

        p_rms_bg = ...
            20e-6 * ...
            10^(Lp_bg(j)/20);

        p_peak_bg = ...
            sqrt(2)*p_rms_bg;

        phaseBG = ...
            2*pi*rand;

        toneBG = ...
            p_peak_bg * ...
            sin( ...
                2*pi*f(j)*t + ...
                phaseBG);

        sigBackground = ...
            sigBackground + ...
            toneBG;

    end

    %% Add small broadband random noise

    broadbandNoise = ...
        2e-5 * randn(size(t));

    sig = ...
        sigSource + ...
        sigBackground + ...
        broadbandNoise;

    %% Save file

    dist = distances(i);
    line = angles(i);

    fileName = ...
        sprintf('mic_%03d.mat',i);

    save( ...
        fullfile(outputFolder,fileName), ...
        'sig', ...
        'fs', ...
        'dist', ...
        'line');

    fprintf( ...
        'Created %s | d = %.1f m | h = %.1f m | angle = %.0f deg\n', ...
        fileName, ...
        dist, ...
        h_r, ...
        line);

end

fprintf('\nSynthetic experiment created successfully.\n')
fprintf('Folder: %s\n',outputFolder)