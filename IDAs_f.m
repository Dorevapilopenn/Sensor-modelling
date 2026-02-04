function [Ceq, f, Cc, a] =IDAs_function(const, Cs, C0)
    spec_names = {'H' 'G' 'D' 'HD' 'HG'}; % species names
    Model      = [ 1   0   0   1    1; ...  % H
                   0   1   0   0    1; ...  % G
                   0   0   1   1    0];     % D
    absorbing  = [ 0   0   1   1    0];
    
    beta = [0, 0, 0, const(1), const(2);
            0, 0, 0, const(1), const(3)];

    beta_f= 10.^beta;
    nsamp= 90; % number of samples
    C_tot_A= Cs*ones(nsamp,1);
    C_tot_B= [C0*(sqrt(.1)*10.^(rand(nsamp,1))), C0*(sqrt(.1)*10.^(rand(nsamp,1)))]; % G concentration
    C_tot_C= Cs*ones(nsamp,1);
    C_eq  = cell(1, 2);
    for i=1:2
        clear C
        C_tot= [C_tot_A C_tot_B(:, i) C_tot_C];   % G concentration
        for j=1:nsamp
            c_comp_guess = [1e-10 1e-10 1e-10 ];
            C(j,:)=NewtonRaphson(Model,beta_f(i, :),C_tot(j,:),c_comp_guess,j);
        end
        C_eq{i}=C;
    end

    % Generation of spectra
    A = spectra_IDA();
    % Simulate the measurements
    randn('state',0);
    D_meas = cell(1,2);
    C_C=cell(1,2);
    for i=1:2
        C_colored= C_eq{i}(:, find(absorbing));
        D_calc = C_colored*A{1};
        D_meas{i} = D_calc + 0.005*max(D_calc(:))*randn(size(D_calc));
        C_C{i} = C_colored;
    end 
    f= D_meas;
    Cc= C_C;
    a= absorbing;
    Ceq= C_eq;
end