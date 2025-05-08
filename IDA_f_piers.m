function [f, Cc, a] =IDA_function(values, c_0, ph)
    spec_names = {'A' 'B' 'C' 'CH' 'AB' 'AC' 'ACH'}; % species names
    Model      = [ 1   0   0   0    1    1    1    ; ...  % A
                   0   1   0   0    1    0    0    ; ...  % B
                   0   0   1   1    0    1    1    ];     % C
    absorbing  = [ 0   0   1   1    0    1    1    ];
    
    beta = [1 1 1  -6.72-ph 5.43 4.91 -2.77-ph;...
            1 1 1  -6.72-ph 3.01 4.91 -2.77-ph;...
            1 1 1  -6.72-ph 4.88 4.65 -3.14-ph ;...
            1 1 1  -6.72-ph 3.03 4.65 -3.14-ph ];

    beta_f= 10.^beta;

    nsamp= 30;
    ncomp= 4;   
    C_tot_A=c_0(1,1)*ones(nsamp,1);
    C_tot_B=[c_0(1,2)*rand(nsamp,1), c_0(1,3)*rand(nsamp,1)]; % G concentration
    C_tot_C=c_0(1,4)*ones(nsamp,1);
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
    lam = 400:1:700;                  
    mean  = [values(7,1) values(8,1) values(9,1) values(10,1); values(7,1) values(8,1) values(11,1) values(12,1)];
    height= [values(13,1) values(14,1) values(15,1) values(16,1); values(13,1) values(14,1) values(17,1) values(18,1)];
    width = [values(19,1) values(20,1) values(21,1) values(22,1); values(19,1) values(20,1) values(23,1) values(24,1)];
    A = cell(1,2);
    for i=1:2
        for j=1:4
            A{i}(j,:)= height(i,j)*gauss(lam,mean(i,j),width(i,j));
        end
    end
    sig_R = 0.000;
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
        D_meas{i} = D_calc + sig_R*randn(size(D_calc));
        C_C{i} = C_colored;
    end 
    f= D_meas;
    Cc= C_C;
    a= absorbing;
end