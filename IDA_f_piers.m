function [Ceq, f, Cc, a, spec] =IDA_function(c_0, ph, ph2)
    spec_names = {'A' 'B' 'C' 'CH' 'AB' 'AC' 'ACH'}; % species names
    Model      = [ 1   0   0   0    1    1    1    ; ...  % A
                   0   1   0   0    1    0    0    ; ...  % B
                   0   0   1   1    0    1    1    ];     % C
    absorbing  = [ 0   0   1   1    0    1    1    ];
    
    beta = [0, 0, 0,  6.72-ph, 8, 6, 10.67-ph;...
            0, 0, 0,  6.72-ph, 2.5, 6, 10.67-ph;...
            0, 0, 0,  6.72-ph2, 2.5, 2.5, 10.3-ph2 ;...
            0, 0, 0,  6.72-ph2, 6, 2.5, 10.3-ph2];

    beta_f= 10.^beta;6
    nsamp= 90; % number of samples
    C_tot_A=[c_0(1,1)*ones(nsamp,1), c_0(1,2)*ones(nsamp,1)];
    C_tot_B=[c_0(1,3)*(sqrt(10)*10.^(rand(nsamp,1))), c_0(1,4)*(sqrt(10)*10.^(rand(nsamp,1)))]; % G concentration
    C_tot_C=c_0(1,5)*ones(nsamp,1);
    C_eq  = cell(1, 4);
    for i=1:4
        clear C
        An= floor((i+1)/2);
        if rem(i,2)==0
            m=2;
        else
            m=1;
        end
        C_tot= [C_tot_A(:, An) C_tot_B(:, m) C_tot_C];   % G concentration
        for j=1:nsamp
            c_comp_guess = [1e-10 1e-10 1e-10 ];
            C(j,:)=NewtonRaphson(Model,beta_f(i, :),C_tot(j,:),c_comp_guess,j);
        end
        C_eq{i}=C;
    end

    % Generation of spectra
    A = piers_spectra_IDA();
    % Simulate the measurements
    randn('state',0);
    D_meas = cell(1,4);
    C_C=cell(1,4);
    for i=1:4
        C_colored= C_eq{i}(:, find(absorbing));
        if i<3
            D_calc = C_colored*A{1};
        else
            D_calc = C_colored*A{2};
        end
        D_meas{i} = D_calc + 0.005*max(D_calc(:))*randn(size(D_calc));
        C_C{i} = C_colored;
    end 
    spec = A;
    f= D_meas;
    Cc= C_C;
    a= absorbing;
    Ceq= C_eq;
end