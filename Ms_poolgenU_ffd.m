function A=Ms_poolgenU_ffd(HD, deltaD, deltaK, deltaC, deltaG)
    % Generate all combinations of parameters for Ms FFD with uncertainty
    % Ensure column vectors
    HD = HD(:);
    deltaD = deltaD(:);
    deltaK = deltaK(:);
    deltaC = deltaC(:);
    deltaG = deltaG(:);

    % Generate all combinations using ndgrid
    [A_grid, B_grid, C_grid, D_grid, E_grid] = ndgrid(HD, deltaD, deltaK, deltaC, deltaG);
    A_vec = A_grid(:);
    B_vec = B_grid(:);
    C_vec = C_grid(:);
    D_vec = D_grid(:);
    E_vec = E_grid(:);



    K1= A_vec;
    K2=K1 + B_vec;
    K3=K2 + C_vec;
    K4=K1 + D_vec;
    K5=K4 + E_vec;

    M=[A_vec, B_vec, C_vec, D_vec, E_vec];

    B_out=[K1, K2, K3, K4, K5, M];
    C_out=[K1, K2, K3, K5, K4, M];
    A=[B_out; C_out];
end 