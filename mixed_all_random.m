tic
clear, clc, close all
spec_names = {'H' 'G' 'D' 'HG1' 'HD'  'DG1'};
Model      = [ 1   0   0    1    1    0; ...  % H
               0   1   0    1    0    1; ...  % G
               0   0   1    0    1    1];     % D


absorbing  = [ 0   0   1   0      1  1];
for p=1:1
    for i=1:4
        log_beta   = [ 0   0   0];
        log_beta(end+1)= 2.75 + 6*rand;
        log_beta(end+1)= 3.2 + 6*rand;
        log_beta(end+1)= 0.5 + 6*rand;
        % log_beta   = [ 0   0   0  5.75  6.2  3.5];
        beta =10.^log_beta;
        beta_f(i, :)=beta;
    end

    nsamp=30;
    c_0     = [16.5e-6 0.0001 16.5e-6];            % tot conc in initial solution Bi,Cl
    ncomp   = length(c_0);   
    C_f = cell(1,2);       % number of components
    for i=1:2
        C_tot_H=c_0(1,1)*ones(nsamp,1);
        C_tot_G=c_0(1,2)*rand(nsamp,1);
        C_tot_D=c_0(1,3)*ones(nsamp,1);
        C_f{i} = [C_tot_H C_tot_G  C_tot_D];
    end
    C_eq  = cell(1, 4);
    for i=1:4
        clear C
        m= floor((i+1)/2);
        C_tot= C_f{m};
        c_comp_guess = [1e-10 1e-10 1e-10 ];  % initial guess for [H],[G] and [D]
        for j=1:nsamp
            C(j,:)=NewtonRaphson(Model,beta_f(i, :),C_tot(j,:),c_comp_guess,j);
            c_comp_guess=C(j,1:ncomp);  % new guess = previous comp conc // why?
        end
        C_eq{i}=C;
    end

                                    %generation of spectra
    lam=400:1:600;                  
    mean  =[450 480 490] ;
    width =[70 60  65];
    height=[50000 70000 95000];
    for i=1:length(mean)
        A(i,:)=height(i)*gauss(lam,mean(i),width(i));
    end
    sig_R=0.000;
    randn('state',0);
    D_meas = cell(1,4);
    for i=1:4
        C_colored= C_eq{i}(:,find(absorbing));
        D_calc = C_colored*A;
        D_meas{i} = D_calc+sig_R*randn(size(D_calc));        %leave only absorbing D, HD and GD species
    end 

    D1= [D_meas{1} D_meas{2}];
    D2= [D_meas{3} D_meas{4}];
    D = [D1; D2];

    [U,S,v]=svd(D,0);
    Sc=U(:,1:3)*S(1:3,1:3);
    %PCAd{p} = Sc
    figure(p)
    plot3(Sc(1:30,1),Sc(1:30,2), Sc(1:30,3),'r*')
    hold on
    plot3(Sc(31:60,1),Sc(31:60,2), Sc(31:60,3),'b*')
end


