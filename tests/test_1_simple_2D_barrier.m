clear
clc
close all

f = [63 125 250 500 1000 2000 4000 8000];
L_W = [95 98 100 102 103 101 97 92];

T_C = 20;
RH = 70;
p = 101325;
c = 343;

d_p = 300;
plateauLength = 20;

receiverDrop = [50 100];

sourceHeight = 2;
receiverHeight = 2;

G_s = 1;
G_m = 1;
G_r = 1;

C0 = 3;

D_source = zeros(size(f));
nReflectingSurfaces = 0;
D_c = calc_Dc(D_source,nReflectingSurfaces);

A_weight = [-26.2 -16.1 -8.6 -3.2 0 1.2 1.0 -1.1];
L_WA = 10*log10(sum(10.^((L_W + A_weight)/10)));

xTerrain = 0:1:d_p;

ResultsAll = struct;

figure
hold on

for k = 1:length(receiverDrop)

    deltaH = receiverDrop(k);

    groundSource = 0;
    groundReceiver = -deltaH;

    zTerrain = zeros(size(xTerrain));

    descending = xTerrain > plateauLength;

    zTerrain(descending) = groundReceiver .* ...
        (xTerrain(descending) - plateauLength) ./ ...
        (d_p - plateauLength);

    S = [0 groundSource + sourceHeight];
    R = [d_p groundReceiver + receiverHeight];

    h_s = sourceHeight;
    h_r = receiverHeight;

    [edges,pathPoints,~,~] = ...
        find_terrain_diffraction_Edges_2D(S,R,xTerrain,zTerrain);

    [z,d_ss,d_sr,e,d,L_diff] = ...
        z_multiEdge_2D(S,R,edges);

    A_div = calc_Atten_div(d);

    [A_atm,alpha] = ...
        calc_Atten_atm(f,T_C,RH,p,d);

    [A_gr,A_s,A_r,A_m] = ...
        calc_Atten_gr(f,h_s,h_r,d_p,G_s,G_m,G_r);

    [A_bar,D_z,K_met,C3,z_min] = ...
        calc_Atten_bar(f,c,z,e,d_ss,d_sr,d,A_gr,'top');

    C_met = calc_Cmet(h_s,h_r,d_p,C0);

    [L_p_DW,L_A_DW,L_A_LT] = ...
        calc_TotalLevel(L_W,D_c,A_div,A_atm,A_gr,A_bar,C_met);

    A_total = A_div + A_atm + A_gr + A_bar;

    A_overall_DW = L_WA - L_A_DW;
    A_overall_LT = L_WA - L_A_LT;

    Results = table(f.',A_atm.',A_gr.',D_z.',A_bar.',A_total.',L_p_DW.', ...
        'VariableNames',{'f_Hz','A_atm_dB','A_gr_dB','D_z_dB','A_bar_dB','A_total_dB','L_p_DW_dB'});

    ResultsAll(k).drop = deltaH;
    ResultsAll(k).S = S;
    ResultsAll(k).R = R;
    ResultsAll(k).edges = edges;
    ResultsAll(k).z = z;
    ResultsAll(k).d = d;
    ResultsAll(k).L_diff = L_diff;
    ResultsAll(k).A_div = A_div;
    ResultsAll(k).A_total = A_total;
    ResultsAll(k).L_A_DW = L_A_DW;
    ResultsAll(k).L_A_LT = L_A_LT;
    ResultsAll(k).A_overall_DW = A_overall_DW;
    ResultsAll(k).A_overall_LT = A_overall_LT;
    ResultsAll(k).A_bar = A_bar;
    ResultsAll(k).D_z = D_z;
    ResultsAll(k).A_gr = A_gr;

    fprintf('\n============================================\n')
    fprintf('RECEIVER %.0f m BELOW SOURCE\n',deltaH)
    fprintf('============================================\n')

    fprintf('S = [%.3f, %.3f] m\n',S(1),S(2))
    fprintf('R = [%.3f, %.3f] m\n',R(1),R(2))

    fprintf('\nGEOMETRY\n')
    fprintf('Horizontal distance = %.3f m\n',d_p)
    fprintf('Direct distance = %.6f m\n',d)
    fprintf('Diffracted path length = %.6f m\n',L_diff)
    fprintf('Path difference z = %.6f m\n',z)

    fprintf('\nDIFFRACTION EDGES\n')
    disp(edges)

    fprintf('A_div = %.6f dB\n',A_div)
    fprintf('C_met = %.6f dB\n',C_met)

    fprintf('\nRESULTS BY OCTAVE BAND\n')
    disp(Results)

    fprintf('L_A_DW = %.6f dBA\n',L_A_DW)
    fprintf('L_A_LT = %.6f dBA\n',L_A_LT)
    fprintf('Overall A-weighted propagation reduction DW = %.6f dB\n',A_overall_DW)
    fprintf('Overall A-weighted propagation reduction LT = %.6f dB\n',A_overall_LT)

    plot(xTerrain,zTerrain,'LineWidth',1.5)
    plot(pathPoints(:,1),pathPoints(:,2),'--','LineWidth',1.5)
    plot(S(1),S(2),'s','MarkerSize',8,'LineWidth',1.5)
    plot(R(1),R(2),'s','MarkerSize',8,'LineWidth',1.5)

end

fprintf('\n============================================\n')
fprintf('COMPARISON\n')
fprintf('============================================\n')

fprintf('50 m drop:\n')
fprintf('z = %.6f m\n',ResultsAll(1).z)
fprintf('A_div = %.6f dB\n',ResultsAll(1).A_div)
fprintf('A_bar = ')
fprintf('%.3f ',ResultsAll(1).A_bar)
fprintf('dB\n')
fprintf('L_A_DW = %.6f dBA\n',ResultsAll(1).L_A_DW)
fprintf('Overall reduction DW = %.6f dB\n',ResultsAll(1).A_overall_DW)

fprintf('\n100 m drop:\n')
fprintf('z = %.6f m\n',ResultsAll(2).z)
fprintf('A_div = %.6f dB\n',ResultsAll(2).A_div)
fprintf('A_bar = ')
fprintf('%.3f ',ResultsAll(2).A_bar)
fprintf('dB\n')
fprintf('L_A_DW = %.6f dBA\n',ResultsAll(2).L_A_DW)
fprintf('Overall reduction DW = %.6f dB\n',ResultsAll(2).A_overall_DW)

fprintf('\nDifference 100 m - 50 m:\n')
fprintf('Delta z = %.6f m\n',ResultsAll(2).z - ResultsAll(1).z)
fprintf('Delta A_div = %.6f dB\n',ResultsAll(2).A_div - ResultsAll(1).A_div)

fprintf('Delta A_bar = ')
fprintf('%.3f ',ResultsAll(2).A_bar - ResultsAll(1).A_bar)
fprintf('dB\n')

fprintf('Delta L_A_DW = %.6f dB\n',ResultsAll(2).L_A_DW - ResultsAll(1).L_A_DW)
fprintf('Additional overall reduction = %.6f dB\n',ResultsAll(2).A_overall_DW - ResultsAll(1).A_overall_DW)

AbarComparison = table( ...
    f.', ...
    ResultsAll(1).A_bar.', ...
    ResultsAll(2).A_bar.', ...
    (ResultsAll(2).A_bar - ResultsAll(1).A_bar).', ...
    'VariableNames',{ ...
        'f_Hz', ...
        'A_bar_50m_dB', ...
        'A_bar_100m_dB', ...
        'Delta_A_bar_dB' ...
    });

fprintf('\nA_bar comparison by octave band:\n')
disp(AbarComparison)