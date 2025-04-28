
clear, clc
load D:/Sens1_D.mat 
D_t=[D_meas1;D_meas2];
[U,S,v]=svd(D_t,0);
Sc=U(:,1:3)*S(1:3,1:3);
figure(3)
plot3(Sc(1:50,1),Sc(1:50,2), Sc(1:50,3),'r*')
hold on
plot3(Sc(51:100,1),Sc(51:100,2), Sc(1:50,3),'b*')
