function A=Ma_poolgenU_ffd(meanHD, deltaD1, deltaK1, deltaC1, deltaG1, deltaHD, deltaD2, deltaK2, deltaC2, deltaG2)
    % Ensure column vectors
    meanHD = meanHD(:);
    deltaD1 = deltaD1(:);
    deltaK1 = deltaK1(:);
    deltaC1 = deltaC1(:);
    deltaG1 = deltaG1(:);
    deltaHD = deltaHD(:);
    deltaD2 = deltaD2(:);
    deltaK2 = deltaK2(:);
    deltaC2 = deltaC2(:);
    deltaG2 = deltaG2(:);

    % Generate all combinations using ndgrid
    [A_grid, B_grid, C_grid, D_grid, E_grid, F_grid, G_grid, H_grid, I_grid, J_grid] = ndgrid(meanHD, deltaD1, deltaK1, deltaC1, deltaG1, deltaHD, deltaD2, deltaK2, deltaC2, deltaG2);

    A_vec = A_grid(:);
    B_vec = B_grid(:);
    C_vec = C_grid(:);
    D_vec = D_grid(:);
    E_vec = E_grid(:);
    F_vec = F_grid(:);
    G_vec = G_grid(:);
    H_vec = H_grid(:);
    I_vec = I_grid(:);
    J_vec = J_grid(:);


    K1= A_vec - F_vec/2;
    K6= A_vec + F_vec/2;
    K2=K1 + B_vec;
    K3=K2 + C_vec;
    K4=K1 + D_vec;
    K5=K4 + E_vec;
    K7=K6 + G_vec;
    K8=K7 + H_vec;
    K9=K6 + I_vec;
    K10=K9 + J_vec;

    M=[A_vec, B_vec, C_vec, D_vec, E_vec, F_vec, G_vec, H_vec, I_vec, J_vec];

    B_out=[K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, M];
    C_out=[K1, K2, K3, K4, K5, K6, K7, K8, K10, K9, M];
    D_out=[K1, K2, K3, K4, K5, K6, K8, K7, K9, K10, M];
    E_out=[K1, K2, K3, K4, K5, K6, K8, K7, K10, K9, M];
    F_out=[K1, K2, K3, K5, K4, K6, K7, K8, K9, K10, M];
    G_out=[K1, K2, K3, K5, K4, K6, K7, K8, K10, K9, M];
    H_out=[K1, K2, K3, K5, K4, K6, K8, K7, K9, K10, M];
    I_out=[K1, K2, K3, K5, K4, K6, K8, K7, K10, K9, M];
    A=[B_out; C_out; D_out; E_out; F_out; G_out; H_out; I_out];
end 