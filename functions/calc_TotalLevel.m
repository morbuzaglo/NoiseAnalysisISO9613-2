% For f = [63,125,250,500,1000,2000,4000,8000] Hz

function [L_p_DW, L_A_DW, L_A_LT] = calc_TotalLevel(L_W, D_c, A_div, A_atm, A_gr, A_bar, C_met)

    A_weight = [-26.2 -16.1 -8.6 -3.2 0.0   1.2   1.0 -1.1];

    L_p_DW = L_W + D_c - A_div - A_atm - A_gr -A_bar;

    L_pA_DW = L_p_DW + A_weight;

    L_A_DW = 10*log10(sum(10.^(L_pA_DW/10)));

    L_A_LT = L_A_DW - C_met;

end