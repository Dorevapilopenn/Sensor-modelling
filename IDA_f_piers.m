function [Ceq, f, Cc, a] =IDA_function(values, c_0, ph, ph2)
    spec_names = {'A' 'B' 'C' 'CH' 'AB' 'AC' 'ACH'}; % species names
    Model      = [ 1   0   0   0    1    1    1    ; ...  % A
                   0   1   0   0    1    0    0    ; ...  % B
                   0   0   1   1    0    1    1    ];     % C
    absorbing  = [ 0   0   1   1    0    1    1    ];
    
    beta = [0, 0, 0,  6.72-ph, 5.46, 2.99, 10.67-ph;...
            0, 0, 0,  6.72-ph, 3.01, 2.99, 10.67-ph;...
            0, 0, 0,  6.72-ph2, 4.88, 2.51, 10.3-ph2 ;...
            0, 0, 0,  6.72-ph2, 3.03, 2.51, 10.3-ph2];

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
    mean  = [499 448 468 451; 499 448 466 452];
    height= [69800 46400 35600 41500; 69800 46400 36100 44300 ];
    width = [27.76 21.49 43.91 22.68; 27.76 21.49 41.66 23.33];
    A = cell(1,2);
    for i=1:2
        for j=1:4
            A{i}(j,:)= height(i,j)*gauss(lam,mean(i,j),width(i,j));
        end
    end
    sig_R = 0.001;
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
    Ceq= C_eq;
end