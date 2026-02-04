function A=IDAs_poolgenU_ffd(HD, deltaD, deltaK)
% Ensure column vectors
HD = HD(:);
deltaD = deltaD(:);
deltaK = deltaK(:);
% Generate all combinations using ndgrid
[A_grid, B_grid, C_grid] = ndgrid(HD, deltaD, deltaK);

A_vec = A_grid(:);
B_vec = B_grid(:);
C_vec = C_grid(:);

K1=A_vec;
K2=K1 + B_vec;
K3=K2 + C_vec;
M=[A_vec, B_vec, C_vec];

A=[K1, K2, K3, M];
end