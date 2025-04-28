
clear, clc
load D:/SensG1.mat 
D_t=[D_meas1 D_meas2];
[U,S,v]=svd(D_t,0);
Sc=U(:,1:3)*S(1:3,1:3);
figure(3)
plot3(Sc(1:50,1),Sc(1:50,2),Sc(1:50,3),'b*')


for i=1:4
    log_beta   = [ 0   0   0]
    log_beta(end+1)= 5.65 + 0.2*rand
    log_beta(end+1)= 6.1 + 0.2*rand
    log_beta(end+1)= 3.4 + 0.2*rand
    % log_beta   = [ 0   0   0  5.75  6.2  3.5];
    beta =10.^-log_beta;
    beta_f(i, :)=beta
end

[coeff, Sc, latent] = pca(D)

scoreA_proj = Sc(1:50, :)
scoreB_proj = Sc(51:100, :)

subplot(3,3,i)
% Scatter plots
scatter(scoreA_proj(:,1), scoreA_proj(:,2), 'bo', 'DisplayName', 'System A');
scatter(scoreB_proj(:,1), scoreB_proj(:,2), 'ro', 'DisplayName', 'System B');

xlabel('PC1');
ylabel('PC2');
title('PCA Projection with Shared Principal Components');
legend;
axis equal;


[U,S,v]=svd(D,0);
    Sc=U(:,1:3)*S(1:3,1:3);
    %PCAd{p} = Sc
    figure(p)
    plot3(Sc(1:50,1),Sc(1:50,2), Sc(1:50,3),'r*')
    hold on
    plot3(Sc(51:100,1),Sc(51:100,2), Sc(51:100,3),'b*')