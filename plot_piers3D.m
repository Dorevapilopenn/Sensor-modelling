% filepath: c:\Users\Propanone\Desktop\Programming\Sensor\Modelling\Sensor-modelling\plot_piers.m
C_0 = [2e-3 1.5e-4 1e-6 1e-5];  % conc A, B1, B2, C respectively
[C_eq, IDA_D, IDA_c, Ia] = IDA_f_piers(values, C_0, 7.2, 7.6);
lgnd= ["C", "CH+", "AC", "ACH+"];

% Plotting

fig1 = figure(1);
% Create a tiled layout (2x3) with compact spacing.
t = tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

nexttile;
plot(IDA_c{1}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("S1G1");

nexttile;
plot(IDA_c{2}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("S1G2");

nexttile;
D = [IDA_D{1}; IDA_D{2}];
[U, S, v] = svd(D, 0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Single(S1)']);

nexttile;
plot(IDA_c{3}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("S2G1");

nexttile;
plot(IDA_c{4}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("S2G2");

nexttile;
D = [IDA_D{3}; IDA_D{4}];
[U, S, v] = svd(D, 0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Single(S2)']);

fig2 = figure(2);
D = [IDA_D{1}, IDA_D{3}; IDA_D{2},IDA_D{4}];
[U, S, v] = svd(D, 0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Array']);