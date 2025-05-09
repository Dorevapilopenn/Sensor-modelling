
clear, clc, close all
spec_names = {'A' 'B' 'C' 'AB' 'AC'}; % species names
Model      = [ 1   0   0   1    1; ...  % A
               0   1   0   1    0; ...  % B
               0   0   1   0    1];     % C
beta =       [ 1   1   1   20   30];
C_tot = [1 .5 .25]; % total concentration in initial solution A,B,C
c_comp_guess = [1e-10 1e-10 1e-10 ];
c = NewtonRaphson(Model, beta, C_tot, c_comp_guess, 1)