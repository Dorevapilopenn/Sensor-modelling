
clear, clc, close all
spec_names = {'H' 'G' 'D' 'GD' 'HG' 'HD' };
Model      = [ 1   0   0   0    1    1   ; ...  % H
               0   1   0   1    1    0; ...  % G
               0   0   1   1    0    1];     % D
absorbing  = [ 0   0   1   1    0    1];
param      = slider_gui_mixed();
names      = param(:,1);
values     = str2double(param(:, 2));
fl=0;
% Generation of spectra
lam   = 400:1:700;                  
mean  = [values(11,1) values(13,1) values(12,1);...
         values(11,1) values(13,1) values(14,1);...
         values(11,1) values(15,1) values(12,1);...
         values(11,1) values(15,1) values(14,1)];
height = [values(16,1) values(18,1) values(17,1);...
          values(16,1) values(18,1) values(19,1);...
          values(16,1) values(20,1) values(17,1);...
          values(16,1) values(20,1) values(19,1)];
width  = [values(21,1) values(23,1) values(22,1);...
          values(21,1) values(23,1) values(24,1);...
          values(21,1) values(25,1) values(22,1);...
          values(21,1) values(25,1) values(24,1)];
A = cell(1,4);
for i = 1:4
    for j = 1:3
        A{i}(j,:) = height(i,j) * gauss(lam, mean(i,j), width(i,j));
    end
end
for p = 1:1000
    beta = [0 0 0 values(3,1) values(2,1) values(1,1);...
            0 0 0 values(8,1) values(7,1) values(6,1);...
            0 0 0 values(5,1) values(4,1) values(1,1);...
            0 0 0 values(10,1) values(9,1) values(6,1)];
    beta_f = 10.^beta;

    nsamp  = 30;
    c_0    = [16.5e-6 0.0001 16.5e-6];  % total conc in initial solution Bi,Cl
    ncomp  = length(c_0);   
    C_f    = cell(1,2); 
    for i = 1:2
        C_tot_H = c_0(1,1)*ones(nsamp,1);
        C_tot_G = c_0(1,2)*rand(nsamp,1);
        C_tot_D = c_0(1,3)*ones(nsamp,1);
        C_f{i}  = [C_tot_H C_tot_G C_tot_D];
    end
    C_eq = cell(1, 4);
    for i = 1:4
        clear C
        m = floor((i+1)/2);
        C_tot = C_f{m};
        c_comp_guess = [1e-10 1e-10 1e-10];
        for j = 1:nsamp
            C(j,:) = NewtonRaphson(Model, beta_f(i, :), C_tot(j,:), c_comp_guess, j); %#ok<SAGROW>
            c_comp_guess = C(j,1:ncomp);
        end
        C_eq{i} = C;
    end

    
    sig_R = 0.000;
    randn('state',0);
    D_meas = cell(1,4);
    for i = 1:4
        C_colored = C_eq{i}(:, find(absorbing));
        D_calc = C_colored * A{i};
        D_meas{i} = D_calc + sig_R * randn(size(D_calc));
    end 

    D1 = [D_meas{1} D_meas{2}];
    D2 = [D_meas{3} D_meas{4}];
    D  = [D1; D2];
    [~, score] = pca(D);
    PC1= score(1:30, :);
    PC2= score(31:60, :);
    f(p, :) = pcaSep(PC1, PC2);
    flp=0;
    if f(p, 1)<3
        flp= flp+1;
    end
    if f(p, 2)<7.8
        flp= flp+1;
    end
    if f(p, 3)>0.01
        flp= flp+1;
    end
    if flp>1
        fl= fl+1;
    end
end
f = transpose(f);
clear mean
for i=1:3
    display(mean(f(i, :)))
    display(std(f(i, :)))
end