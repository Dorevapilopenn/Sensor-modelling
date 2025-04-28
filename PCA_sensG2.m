
clear, clc
load D:/SensG2.mat 
D_t=[D_meas1 D_meas2];
[U,S,v]=svd(D_t,0);
Sc=U(:,1:3)*S(1:3,1:3);
figure(3)
plot3(Sc(1:50,1),Sc(1:50,2),Sc(1:50,3),'b*')