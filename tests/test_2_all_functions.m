clear
clc
close all

%% ============================================================
%  PATHS
% =============================================================

addpath(fullfile(pwd,'functions'));

experimentFolder = fullfile(pwd,'experiment');

%% ============================================================
%  GENERAL PARAMETERS
% =============================================================

f = [63 125 250 500 1000 2000 4000 8000];

T_C = 20;
RH = 70;
p = 101325;
c = 343;

h_s = 1.5;

%% ============================================================
%  EXPERIMENT MICROPHONE GEOMETRY
% =============================================================

% Mic:
% 1   2   3   4   5    6   7   8   9   10  11  12  13

d_p_all = ...
    [7 14 21 50 250 7 14 21 50 7 14 21 50];

h_r_all = ...
    [1.5 1.5 1.5 1.5 -11 ...
     1.5 1.5 1.5 1.5 ...
     1.5 1.5 1.5 1.5];

nMics = length(d_p_all);
nFreq = length(f);

% Fallback angles if "line" does not exist in the MAT files
defaultLineAnglesDeg = [0 120 240];

theta_all = [
    repmat(defaultLineAnglesDeg(1),1,5) ...
    repmat(defaultLineAnglesDeg(2),1,4) ...
    repmat(defaultLineAnglesDeg(3),1,4)
];

%% Ground properties

G_s_all = ones(1,nMics);
G_m_all = ones(1,nMics);
G_r_all = ones(1,nMics);

%% Reflecting-surface correction

D_surface_all = zeros(nMics,nFreq);

%% ============================================================
%  ANALYSIS TIME WINDOWS
% =============================================================

fprintf('\n============================================\n')
fprintf('SELECT ANALYSIS WINDOWS\n')
fprintf('============================================\n')

t1 = input('Source activity start time [s]: ');
t2 = input('Source activity end time [s]: ');

t1_bg = input('Background start time [s]: ');
t2_bg = input('Background end time [s]: ');

%% ============================================================
%  LOAD ALL 13 MICROPHONES
% =============================================================

Lp_total_all = nan(nMics,nFreq);
Lp_bg_all = nan(nMics,nFreq);

fs_all = nan(nMics,1);

fprintf('\n============================================\n')
fprintf('LOAD EXPERIMENT DATA\n')
fprintf('============================================\n')

for i = 1:nMics

    fileName = sprintf('mic_%03d.mat',i);
    filePath = fullfile(experimentFolder,fileName);

    if ~isfile(filePath)
        error('Missing file: %s',filePath);
    end

    data = load(filePath);

    if ~isfield(data,'sig')
        error('%s does not contain sig.',fileName);
    end

    if ~isfield(data,'fs')
        error('%s does not contain fs.',fileName);
    end

    sig = data.sig(:);
    fs = data.fs;

    fs_all(i) = fs;

    %% Use experiment line angle if available

    if isfield(data,'line')
        theta_all(i) = data.line;
    end

    %% Optional consistency check for distance

    if isfield(data,'dist')

        if abs(data.dist-d_p_all(i)) > 0.5
            warning( ...
                'Mic %d: file dist = %.2f m, geometry dist = %.2f m.', ...
                i,data.dist,d_p_all(i));
        end

    end

    %% Sampling-rate test

    G = 10^(3/10);
    f_high_max = max(f)*sqrt(G);

    if fs <= 2*f_high_max

        error( ...
            'Mic %d: fs = %.0f Hz is too low for the 8 kHz octave band.', ...
            i,fs);

    elseif fs < 2.5*f_high_max

        warning( ...
            'Mic %d: fs = %.0f Hz has limited margin for 8 kHz.', ...
            i,fs);

    end

    %% Calculate octave-band levels

    Lp_total_all(i,:) = ...
        calc_Lp_total( ...
            sig,fs,f,t1,t2,'octave');

    Lp_bg_all(i,:) = ...
        calc_Lp_bg( ...
            sig,fs,f,t1_bg,t2_bg,'octave');

    fprintf( ...
        'Mic %02d | d_p = %6.1f m | h = %6.1f m | theta = %7.1f deg | fs = %.0f Hz\n', ...
        i,d_p_all(i),h_r_all(i),theta_all(i),fs);

end

%% ============================================================
%  TEST calc_Lp_band AND OCTAVE FILTER OPTIONS
% =============================================================

fprintf('\n============================================\n')
fprintf('OCTAVE FILTER TEST\n')
fprintf('============================================\n')

data = load(fullfile(experimentFolder,'mic_001.mat'));

sig = data.sig(:);
fs = data.fs;

fc_test = 1000;

[Lp_test,p_rms_test,sigBand_test] = ...
    calc_Lp_band( ...
        sig,fs,fc_test,t1,t2,'octave');

fprintf('calc_Lp_band at %d Hz:\n',fc_test)
fprintf('Lp = %.3f dB\n',Lp_test)
fprintf('p_rms = %.6e Pa\n',p_rms_test)

i1 = max(1,round(t1*fs)+1);
i2 = min(length(sig),round(t2*fs));

sigSegment = sig(i1:i2);

%% Automatic filter

[sigAuto,methodUsed] = ...
    octaveFilter_auto( ...
        sigSegment,fs,fc_test);

fprintf('octaveFilter_auto method: %s\n',methodUsed)

%% Custom IEC filter

sigIEC = ...
    octaveFilter_IEC( ...
        sigSegment,fs,fc_test);

fprintf('octaveFilter_IEC tested.\n')

%% Signal Processing Toolbox version

if exist('butter','file') == 2

    sigSignal = ...
        octaveFilter_signalToolbox( ...
            sigSegment,fs,fc_test);

    fprintf('octaveFilter_signalToolbox tested.\n')

else

    fprintf('Signal Processing Toolbox option unavailable.\n')

end

%% Audio Toolbox version

if exist('octaveFilter','class') == 8 || ...
   exist('octaveFilter','file') == 2

    sigAudio = ...
        octaveFilter_audioToolbox( ...
            sigSegment,fs,fc_test);

    fprintf('octaveFilter_audioToolbox tested.\n')

else

    fprintf('Audio Toolbox option unavailable.\n')

end

%% Custom filter Class-1 check

[class1OK,fTest,attenuationTest] = ...
    check_octaveFilter_class1( ...
        fs,fc_test);

fprintf('check_octaveFilter_class1 result = %d\n',class1OK)

%% ============================================================
%  SPLIT CALIBRATION AND VALIDATION MICROPHONES
% =============================================================

% Close mics:
% 7, 14, 21 m in all three directions

closeIdx = find(d_p_all <= 21);

% Far mics:
% 50 m in three directions + Mic 5 at 250 m

farIdx = find(d_p_all > 21);

fprintf('\nCalibration mics:\n')
disp(closeIdx)

fprintf('Validation mics:\n')
disp(farIdx)

%% ============================================================
%  CLOSE-MIC GEOMETRY
% =============================================================

d_p_close = d_p_all(closeIdx);
h_r_close = h_r_all(closeIdx);

d_close = sqrt( ...
    d_p_close.^2 + ...
    (h_s-h_r_close).^2);

theta_close = theta_all(closeIdx);

G_s_close = G_s_all(closeIdx);
G_m_close = G_m_all(closeIdx);
G_r_close = G_r_all(closeIdx);

D_surface_close = ...
    D_surface_all(closeIdx,:);

Lp_total_close = ...
    Lp_total_all(closeIdx,:);

Lp_bg_close = ...
    Lp_bg_all(closeIdx,:);

%% ============================================================
%  ESTIMATE LW + SOURCE DIRECTIVITY
% =============================================================

[L_W,L_W_std,L_W_mics,D_source_raw, ...
 Lp_source_close,validBG_close] = ...
    estimate_Lw_from_close_mics( ...
        f, ...
        Lp_total_close, ...
        Lp_bg_close, ...
        d_close, ...
        d_p_close, ...
        h_s, ...
        h_r_close, ...
        G_s_close, ...
        G_m_close, ...
        G_r_close, ...
        T_C, ...
        RH, ...
        p, ...
        D_surface_close, ...
        theta_close);

%% ============================================================
%  AVERAGE DIRECTIVITY OF MICS ON SAME LINE
% =============================================================

calibAnglesDeg = unique(theta_close(:));

D_source = nan( ...
    length(calibAnglesDeg), ...
    nFreq);

for i = 1:length(calibAnglesDeg)

    idx = theta_close(:) == calibAnglesDeg(i);

    D_source(i,:) = ...
        mean( ...
            D_source_raw(idx,:), ...
            1,'omitnan');

end

%% ============================================================
%  DISPLAY SOURCE CHARACTERIZATION
% =============================================================

fprintf('\n============================================\n')
fprintf('ESTIMATED SOURCE CHARACTERIZATION\n')
fprintf('============================================\n')

SourceTable = table( ...
    f.', ...
    L_W.', ...
    L_W_std.', ...
    'VariableNames',{ ...
        'f_Hz', ...
        'L_W_dB', ...
        'Std_dB'});

disp(SourceTable)

DirectivityTable = ...
    array2table( ...
        D_source, ...
        'VariableNames',compose('Hz_%d',f));

DirectivityTable.Angle_deg = ...
    calibAnglesDeg;

DirectivityTable = ...
    movevars( ...
        DirectivityTable, ...
        'Angle_deg', ...
        'Before',1);

disp(DirectivityTable)

%% ============================================================
%  TEST DIRECTIVITY INTERPOLATION
% =============================================================

thetaTest = 45;

D_source_test = ...
    estimate_Dsource_at_angle( ...
        calibAnglesDeg, ...
        D_source, ...
        thetaTest);

fprintf('\nD_source interpolated at %.1f deg:\n',thetaTest)
disp(D_source_test)

%% ============================================================
%  PRECOMPUTE ATTENUATIONS FOR VALIDATION MICS
% =============================================================

nFar = length(farIdx);

d_p_far = d_p_all(farIdx);
h_r_far = h_r_all(farIdx);
theta_far = theta_all(farIdx);

d_far = sqrt( ...
    d_p_far.^2 + ...
    (h_s-h_r_far).^2);

G_s_far = G_s_all(farIdx);
G_m_far = G_m_all(farIdx);
G_r_far = G_r_all(farIdx);

D_surface_far = ...
    D_surface_all(farIdx,:);

Lp_total_far = ...
    Lp_total_all(farIdx,:);

Lp_bg_far = ...
    Lp_bg_all(farIdx,:);

A_div_far = nan(nFar,nFreq);
A_atm_far = nan(nFar,nFreq);
A_gr_far = nan(nFar,nFreq);

%% IMPORTANT:
% Until the real terrain profile of the experiment is supplied,
% validation assumes no barrier/topographic diffraction.

A_bar_far = zeros(nFar,nFreq);

for i = 1:nFar

    A_div_i = ...
        calc_Atten_div( ...
            d_far(i));

    A_div_far(i,:) = ...
        A_div_i .* ones(size(f));

    [A_atm_far(i,:),~] = ...
        calc_Atten_atm( ...
            f,T_C,RH,p,d_far(i));

    [A_gr_far(i,:),~,~,~] = ...
        calc_Atten_gr( ...
            f, ...
            h_s, ...
            h_r_far(i), ...
            d_p_far(i), ...
            G_s_far(i), ...
            G_m_far(i), ...
            G_r_far(i));

end

%% ============================================================
%  VALIDATE FAR MICS
% =============================================================

results = ...
    validate_far_mics( ...
        f, ...
        L_W, ...
        calibAnglesDeg, ...
        D_source, ...
        Lp_total_far, ...
        Lp_bg_far, ...
        theta_far, ...
        D_surface_far, ...
        A_div_far, ...
        A_atm_far, ...
        A_gr_far, ...
        A_bar_far);

%% ============================================================
%  VALIDATION SUMMARY
% =============================================================

fprintf('\n============================================\n')
fprintf('VALIDATION SUMMARY\n')
fprintf('============================================\n')

ValidationTable = table( ...
    farIdx(:), ...
    d_p_far(:), ...
    h_r_far(:), ...
    theta_far(:), ...
    results.Bias, ...
    results.MAE, ...
    results.RMSE, ...
    results.LpA_measured, ...
    results.LpA_predicted, ...
    results.Delta_LA, ...
    'VariableNames',{ ...
        'Mic', ...
        'HorizontalDistance_m', ...
        'ReceiverHeight_m', ...
        'Angle_deg', ...
        'Bias_dB', ...
        'MAE_dB', ...
        'RMSE_dB', ...
        'Measured_dBA', ...
        'Predicted_dBA', ...
        'Delta_dBA'});

disp(ValidationTable)

%% Error by frequency

FrequencyValidationTable = table( ...
    f.', ...
    results.Bias_freq.', ...
    results.MAE_freq.', ...
    results.RMSE_freq.', ...
    'VariableNames',{ ...
        'f_Hz', ...
        'Bias_dB', ...
        'MAE_dB', ...
        'RMSE_dB'});

disp(FrequencyValidationTable)

%% ============================================================
%  VALIDATION PLOT
% =============================================================

colors = lines(nFar);

figure
hold on

for i = 1:nFar

    semilogx( ...
        f, ...
        results.Lp_measured(i,:), ...
        'o-', ...
        'Color',colors(i,:), ...
        'LineWidth',1.5, ...
        'DisplayName', ...
        sprintf('Mic %d measured',farIdx(i)));

    semilogx( ...
        f, ...
        results.Lp_predicted(i,:), ...
        '--', ...
        'Color',colors(i,:), ...
        'LineWidth',1.5, ...
        'DisplayName', ...
        sprintf('Mic %d predicted',farIdx(i)));

end

grid on
xlabel('Frequency [Hz]')
ylabel('L_p [dB]')
title('Validation - measured vs predicted')
legend('Location','best')

%% Error plot

figure
hold on

for i = 1:nFar

    semilogx( ...
        f, ...
        results.Delta_L(i,:), ...
        'o-', ...
        'Color',colors(i,:), ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Mic %d',farIdx(i)));

end

yline(0,'--')

grid on
xlabel('Frequency [Hz]')
ylabel('\DeltaL = measured - predicted [dB]')
title('Validation error')
legend('Location','best')

%% ============================================================
%  BACKGROUND-NOISE DATABASE TEST
% =============================================================

fprintf('\n============================================\n')
fprintf('BACKGROUND-NOISE OPTIONS\n')
fprintf('============================================\n')

%% Typical values

[f_bg_typical,Lp_bg_typical,LpA_bg_typical,descTypical] = ...
    get_background_spectrum_typical('rural');

fprintf('\nTypical background:\n')
fprintf('%s\n',descTypical)
fprintf('LpA = %.2f dBA\n',LpA_bg_typical)

%% Illinois regulatory reference

[f_bg_IL,Lp_bg_IL,LpA_bg_IL,descIL] = ...
    get_background_spectrum_Illinois_regulation( ...
        5,'night');

fprintf('\nIllinois regulation reference:\n')
fprintf('%s\n',descIL)
fprintf('LpA = %.2f dBA\n',LpA_bg_IL)

%% ============================================================
%  NEW TERRAIN PREDICTION
% =============================================================

fprintf('\n============================================\n')
fprintf('NEW TERRAIN SCENARIO\n')
fprintf('============================================\n')

%% Interesting synthetic terrain

xTerrain = linspace(0,500,1001);

zTerrain = ...
      5*exp(-((xTerrain-90)/40).^2) ...
    +14*exp(-((xTerrain-220)/32).^2) ...
    - 6*exp(-((xTerrain-320)/55).^2) ...
    + 9*exp(-((xTerrain-410)/38).^2);

sourceGround = ...
    interp1(xTerrain,zTerrain,0);

receiverGround = ...
    interp1(xTerrain,zTerrain,500);

newReceiverHeight = 1.5;

S = ...
    [0 sourceGround+h_s];

R = ...
    [500 receiverGround+newReceiverHeight];

%% ============================================================
%  FIND 2D TERRAIN DIFFRACTION EDGES
% =============================================================

edges = ...
    find_terrain_diffraction_Edges_2D( ...
        S,R,xTerrain,zTerrain);

fprintf('\nDetected diffraction edges:\n')
disp(edges)

%% Diffraction path geometry

[z,d_ss,d_sr,e,d_new,L_diff] = ...
    z_multiEdge_2D( ...
        S,R,edges);

fprintf('Direct distance  = %.3f m\n',d_new)
fprintf('Diffracted path  = %.3f m\n',L_diff)
fprintf('Path difference  = %.6f m\n',z)

%% ============================================================
%  NEW-SCENARIO ATTENUATIONS
% =============================================================

d_p_new = abs(R(1)-S(1));

A_div_new_scalar = ...
    calc_Atten_div(d_new);

A_div_new = ...
    A_div_new_scalar .* ones(size(f));

[A_atm_new,alpha_new] = ...
    calc_Atten_atm( ...
        f,T_C,RH,p,d_new);

G_s_new = 1;
G_m_new = 1;
G_r_new = 1;

[A_gr_new,~,~,~] = ...
    calc_Atten_gr( ...
        f, ...
        h_s, ...
        newReceiverHeight, ...
        d_p_new, ...
        G_s_new, ...
        G_m_new, ...
        G_r_new);

%% ============================================================
%  DIFFRACTION / BARRIER ATTENUATION
% =============================================================

diffractionType = 'top';

[A_bar_new,D_z_new,K_met_new,C3_new,z_min_new] = ...
    calc_Atten_bar( ...
        f, ...
        c, ...
        z, ...
        e, ...
        d_ss, ...
        d_sr, ...
        d_new, ...
        A_gr_new, ...
        diffractionType);

%% ============================================================
%  NEW RECEIVER DIRECTIVITY
% =============================================================

theta_new = 70;

D_source_new = ...
    estimate_Dsource_at_angle( ...
        calibAnglesDeg, ...
        D_source, ...
        theta_new);

D_surface_new = zeros(size(f));

%% ============================================================
%  PREDICT SOURCE LEVEL AT NEW RECEIVER
% =============================================================

[Lp_source_new,D_source_check] = ...
    predict_Lp_at_mic( ...
        f, ...
        L_W, ...
        theta_new, ...
        calibAnglesDeg, ...
        D_source, ...
        D_surface_new, ...
        A_div_new, ...
        A_atm_new, ...
        A_gr_new, ...
        A_bar_new);

%% ============================================================
%  SELECT BACKGROUND MODEL FOR NEW SCENARIO
% =============================================================

% ------------------------------------------------------------
% OPTION 1 - typical reference spectrum
% ------------------------------------------------------------

backgroundModel = 'typical';

switch lower(backgroundModel)

    case 'typical'

        [~,Lp_bg_new,LpA_bg_reference,bgDescription] = ...
            get_background_spectrum_typical( ...
                'rural');

    case 'illinois'

        [~,Lp_bg_new,LpA_bg_reference,bgDescription] = ...
            get_background_spectrum_Illinois_regulation( ...
                5,'night');

    otherwise

        error( ...
            'backgroundModel must be ''typical'' or ''illinois''.');

end

fprintf('\nBackground model:\n')
fprintf('%s\n',bgDescription)
fprintf('Background = %.2f dBA\n',LpA_bg_reference)

%% ============================================================
%  ADD SOURCE + BACKGROUND ENERGETICALLY
% =============================================================

Lp_total_new = ...
    10*log10( ...
        10.^(Lp_source_new/10) + ...
        10.^(Lp_bg_new/10));

%% A-weighting

A_weight = ...
    [-26.2 -16.1 -8.6 -3.2 0 1.2 1.0 -1.1];

LpA_source_new = ...
    10*log10(sum( ...
        10.^((Lp_source_new+A_weight)/10)));

LpA_background_new = ...
    10*log10(sum( ...
        10.^((Lp_bg_new+A_weight)/10)));

LpA_total_new = ...
    10*log10(sum( ...
        10.^((Lp_total_new+A_weight)/10)));

%% ============================================================
%  NEW-SCENARIO RESULTS
% =============================================================

PredictionTable = table( ...
    f.', ...
    L_W.', ...
    D_source_new.', ...
    A_div_new.', ...
    A_atm_new.', ...
    A_gr_new.', ...
    D_z_new.', ...
    A_bar_new.', ...
    Lp_source_new.', ...
    Lp_bg_new.', ...
    Lp_total_new.', ...
    'VariableNames',{ ...
        'f_Hz', ...
        'L_W', ...
        'D_source', ...
        'A_div', ...
        'A_atm', ...
        'A_gr', ...
        'D_z', ...
        'A_bar', ...
        'Lp_source', ...
        'Lp_background', ...
        'Lp_total'});

disp(PredictionTable)

fprintf('\nSource only   = %.2f dBA\n',LpA_source_new)
fprintf('Background    = %.2f dBA\n',LpA_background_new)
fprintf('Combined      = %.2f dBA\n',LpA_total_new)

%% ============================================================
%  TERRAIN PLOT
% =============================================================

figure

plot(xTerrain,zTerrain,'LineWidth',1.5)
hold on

plot(S(1),S(2),'o','MarkerSize',8)
plot(R(1),R(2),'s','MarkerSize',8)

if ~isempty(edges)

    plot( ...
        edges(:,1), ...
        edges(:,2), ...
        'x', ...
        'MarkerSize',10, ...
        'LineWidth',2);

    pathPoints = ...
        [S; edges; R];

    plot( ...
        pathPoints(:,1), ...
        pathPoints(:,2), ...
        '--', ...
        'LineWidth',1.5);

end

grid on
xlabel('x [m]')
ylabel('Elevation [m]')
title('New prediction terrain')
legend( ...
    'Terrain', ...
    'Source', ...
    'Receiver', ...
    'Diffraction edges', ...
    'Diffracted path', ...
    'Location','best')

%% ============================================================
%  MASS-LAW TEST
% =============================================================

rhoConcrete = 2400;
tConcrete = 0.15;

R_masslaw = ...
    calc_Atten_masslaw( ...
        f, ...
        rhoConcrete, ...
        tConcrete, ...
        'diffuse');

fprintf('\n150 mm concrete mass-law attenuation:\n')
disp(R_masslaw)

%% ============================================================
%  OPTIONAL FUNCTIONS
% =============================================================

% These functions exist in your functions folder but are not necessarily
% required in the present instantaneous experiment workflow:
%
% calc_Cmet
% calc_Dc
% calc_Dsource
% calc_TotalLevel
% combine_diffraction_paths
% find_terrain_diffraction_paths_3D
%
% Reasons:
%
% calc_Cmet:
% long-term meteorological correction, not normally applied to one
% specific measured experiment.
%
% calc_Dc:
% D_surface is supplied explicitly here.
%
% calc_Dsource:
% source directivity is estimated directly from the close microphones.
%
% calc_TotalLevel:
% the equations are written explicitly here so every term can be checked.
%
% combine_diffraction_paths:
% needed when top + left + right diffraction paths are available.
%
% find_terrain_diffraction_paths_3D:
% use when a full X-Y-Z terrain map is available. The example above is 2D.