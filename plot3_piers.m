% filepath: c:\Users\Propanone\Desktop\Programming\Sensor\Modelling\Sensor-modelling\plot_piers.m
C_0 = [2e-3 1e-8 1e-4 1e-5];  % conc A, B1, B2, C respectively
[C_eq, IDA_D, IDA_c, Ia,  A] = IDA_f_piers(values, C_0, 7.2, 7.6);
lgnd= ["C", "CH+", "AC", "ACH+"];

% Plotting

fig1 = figure(1);
% Create a tiled layout (2x3) with compact spacing.
t = tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');


nexttile;
D = [IDA_D{1}; IDA_D{2}];
[U, S, v] = svd(D, 0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Single(S1)']);

nexttile;
D = [IDA_D{3}; IDA_D{4}];
[U, S, v] = svd(D, 0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Single(S2)']);

nexttile;
D = [IDA_D{1}, IDA_D{3}; IDA_D{2},IDA_D{4}];
[U, S, v] = svd(D, 0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Array']);