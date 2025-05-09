function [f, Cc, a] =IDA_function(values, c_0)
    spec_names = {'H' 'G' 'D' 'HG' 'HD'};
    Model      = [ 1   0   0    1    1  ; ...  % H
                   0   1   0    1    0  ; ...  % G
                   0   0   1    0    1  ];     % D
    absorbing  = [ 0   0   1    0    1  ];
    
    beta = [0 0 0 values(2,1) values(1,1);...
            0 0 0 values(3,1) values(1,1);...
            0 0 0 values(5,1) values(4,1);...
            0 0 0 values(6,1) values(4,1)];
    beta_f= 10.^beta;

    nsamp=30;
    ncomp   = 3;   
    C_tot_H=c_0(1,1)*ones(nsamp,1);
    C_tot_G=[c_0(1,2)*rand(nsamp,1), c_0(1,3)*rand(nsamp,1)]; % G concentration
    C_tot_D=c_0(1,4)*ones(nsamp,1);
    C_eq  = cell(1, 4);
    for i=1:4
        clear C
        if rem(i,2)==0
            m=2;
        else
            m=1;
        end
        C_tot= [C_tot_H C_tot_G(:, m) C_tot_D];   % G concentration
        c_comp_guess = [1e-10 1e-10 1e-10 ];
        for j=1:nsamp
            C(j,:)=NewtonRaphson(Model,beta_f(i, :),C_tot(j,:),c_comp_guess,j);
            c_comp_guess=C(j,1:ncomp);
        end
        C_eq{i}=C;
    end

    % Generation of spectra
    lam = 400:1:700;                  
    mean  = [values(7,1) values(8,1); values(7,1) values(9,1)];
    width = [values(13,1) values(14,1); values(13,1) values(15,1)];
    height= [values(10,1) values(11,1); values(10,1) values(12,1)];
    A = cell(1,2);
    for i=1:2
        for j=1:2
            A{i}(j,:)= height(i,j)*gauss(lam,mean(i,j),width(i,j));
        end
    end
    sig_R = 0.0001;
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