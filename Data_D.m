clear, clc, close all
spec_names = {'H' 'G1' 'G2' 'D' 'HG1' 'HG2' 'HD' 'DG1' 'DG2' };
Model      = [ 1   0    0    0   1     1     1    0     0     ; ...  % H
               0   1    0    0   1     0     0    1     0     ; ...  % G1
               0   0    1    0   0     1     0    0     1     ; ...  % G2
               0   0    0    1   0     0     1    1     1];     % D
log_beta   = [ 0   0    0    0  6.25   5.8  5.2   4    4.4];
% log_beta   = [ 0   0   0  5.75  6.2  3.5];

absorbing  = [ 0   0    0    1   0     0     1    1     1];
beta =10.^log_beta;

nsamp=50;
c_0     = [16.5e-6 0.0001 0.0001 16.5e-6];            % tot conc in initial solution Bi,Cl
ncomp   = length(c_0);          % number of components

C_tot_H=c_0(1,1)*ones(nsamp,1);
C_tot_G1=c_0(1,2)*rand(nsamp,1);
C_tot_G2=c_0(1,3)*rand(nsamp,1);
C_tot_D=c_0(1,4)*ones(nsamp,1);
C_tot=[C_tot_H C_tot_G1 C_tot_G2  C_tot_D];

c_comp_guess = [1e-10 1e-10 1e-10 1e-10 ];  % initial guess for [H],[G] and [D]
for i=1:nsamp
    C(i,:)=NewtonRaphson(Model,beta,C_tot(i,:),c_comp_guess,i);
    c_comp_guess=C(i,1:ncomp);  % new guess = previous comp conc // why?
end
C_colored=C(:,find(absorbing));         %leave only absorbing D, HD and GD species

                                   %generation of spectra
lam=400:1:600;                  
mean  =[450 480 490 510] ;
width =[70 60  65 65];
height=[50000 70000 95000 90000];
for i=1:length(mean)
    A(i,:)=height(i)*gauss(lam,mean(i),width(i));
end

D_calc=C_colored*A;
sig_R=0.000;
randn('state',0);
D_meas1=D_calc+sig_R*randn(size(D_calc));

figure(1);
subplot(2,2,1);
plot(lam,A,'.-');
xlabel('Wavelength');ylabel('Molar Absorptivity')

subplot(2,2,3)
plot(C_colored,'.-');
legend(spec_names(find(absorbing)));
xlabel('# Samples');ylabel('[species]')

figure(2);
subplot(1,2,1)
plot(lam,D_meas1,'y');xlabel('Wavelength');ylabel('Absorbance')

%%
Model      = [ 1   0    0    0   1     1     1    0     0     ; ...  % H
               0   1    0    0   1     0     0    1     0     ; ...  % G1
               0   0    1    0   0     1     0    0     1     ; ...  % G2
               0   0    0    1   0     0     1    1     1];     % D
log_beta   = [ 0   0    0    0  5.1   5.2  4.5    5   4.9];
% log_beta   = [ 0   0   0  5.75  6.2  3.5];

absorbing  = [ 0   0    0    1   0     0     1    1     1];
beta =10.^log_beta;

nsamp=50;
c_0     = [16.5e-6 0.0001 0.0001 16.5e-6];            % tot conc in initial solution Bi,Cl
ncomp   = length(c_0);          % number of components

C_tot_H=c_0(1,1)*ones(nsamp,1);
C_tot_G1=c_0(1,2)*rand(nsamp,1);
C_tot_G2=c_0(1,3)*rand(nsamp,1);
C_tot_D=c_0(1,4)*ones(nsamp,1);
C_tot=[C_tot_H C_tot_G1 C_tot_G2  C_tot_D];

c_comp_guess = [1e-10 1e-10 1e-10 1e-10 ];  % initial guess for [H],[G] and [D]
for i=1:nsamp
    C(i,:)=NewtonRaphson(Model,beta,C_tot(i,:),c_comp_guess,i);
    c_comp_guess=C(i,1:ncomp);  % new guess = previous comp conc // why?
end
C_colored=C(:,find(absorbing));         %leave only absorbing D, HD and GD species

                                   %generation of spectra
lam=400:1:600;                  
mean  =[450 480 490 510] ;
width =[70 60  65 65];
height=[50000 70000 95000 90000];
for i=1:length(mean)
    A(i,:)=height(i)*gauss(lam,mean(i),width(i));
end

D_calc=C_colored*A;
sig_R=0.000;
randn('state',0);
D_meas2=D_calc+sig_R*randn(size(D_calc));

figure(1);
subplot(2,2,2);
plot(lam,A,'.-');
xlabel('Wavelength');ylabel('Molar Absorptivity')

subplot(2,2,4)
plot(C_colored,'.-');
legend(spec_names(find(absorbing)));
xlabel('# Samples');ylabel('[species]')

figure(2);
subplot(1,2,2)
plot(lam,D_meas2,'y');xlabel('Wavelength');ylabel('Absorbance')
save D:/Sens1_D.mat D_meas1 D_meas2