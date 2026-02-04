function [Ceq, f, Cc, a] =Ma_function(const, Cs, C0)
    spec_names = {'H' 'G' 'D' 'HD' 'HG' 'GD'}; % species names
    Model =      [ 1   0   0   1    1    0; ...  % H
                   0   1   0   0    1    1; ...  % G
                   0   0   1   1    0    1];     % D
    absorbing =  [ 0   0   1   1    0    1];

    beta = [0, 0, 0, const(1), const(2), const(4);
            0, 0, 0, const(1), const(3), const(5);
            0, 0, 0, const(6), const(7), const(9);
            0, 0, 0, const(6), const(8), const(10)];

    beta_f= 10.^beta;
    nsamp= 90; % number of samples
    C_tot_A= Cs*ones(nsamp,1);
    C_tot_B= [C0*(sqrt(.1)*10.^(rand(nsamp,1))), C0*(sqrt(.1)*10.^(rand(nsamp,1)))]; % G concentration
    C_tot_C= Cs*ones(nsamp,1);
    C_eq  = cell(1, 4);
    for i=1:4
        clear C
        if rem(i,2)==0
            m=2;
        else
            m=1;
        end
        C_tot= [C_tot_A C_tot_B(:, m) C_tot_C];   % G concentration
        for j=1:nsamp
            c_comp_guess = [1e-10 1e-10 1e-10 ];
            C(j,:)=NewtonRaphson(Model,beta_f(i, :),C_tot(j,:),c_comp_guess,j);
        end
        C_eq{i}=C;
    end

    % Generation of spectra
    A = spectra_M();
    % Simulate the measurements
    randn('state',0);
    D_meas = cell(1,4);
    C_C=cell(1,4);
    for i=1:4
        C_colored= C_eq{i}(:, find(absorbing));
        D_calc = C_colored*A{i};
        D_meas{i} = D_calc + 0.005*max(D_calc(:))*randn(size(D_calc));
        C_C{i} = C_colored;
    end
    f= D_meas;
    Cc= C_C;
    a= absorbing;
    Ceq= C_eq;
end