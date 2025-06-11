% filepath: c:\Users\Propanone\Desktop\Programming\Sensor\Modelling\Sensor-modelling\plot_piers.m
param = box_gui_IDA_piers();
names = param(:, 1);
values = str2double(param(:, 2));
C_0 = [16.5e-2 1.2e-5 1.3e-5 16.5e-2];  % conc A, B1, B2, C respectively
[IDA_D, IDA_c, Ia] = IDA_f_piers(values, C_0, 7.2, 7.6);
lgnd= ["C", "CH+", "AC", "ACH+"]

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
Sc = U(:, 1:2) * S(1:2, 1:2);
scatter(Sc(1:30, 1), Sc(1:30, 2), 36, 'r', 'filled');
hold on;
scatter(Sc(31:60, 1), Sc(31:60, 2), 36, 'b', 'filled');
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
Sc = U(:, 1:2) * S(1:2, 1:2);
scatter(Sc(1:30, 1), Sc(1:30, 2), 36, 'r', 'filled');
hold on;
scatter(Sc(31:60, 1), Sc(31:60, 2), 36, 'b', 'filled');
title(['Single(S2)']);

fig2 = figure(3);
D = [IDA_D{1}, IDA_D{3}; IDA_D{2},IDA_D{4}];
[U, S, v] = svd(D, 0);
Sc = U(:, 1:2) * S(1:2, 1:2);
scatter(Sc(1:30, 1), Sc(1:30, 2), 36, 'r', 'filled');
hold on;
scatter(Sc(31:60, 1), Sc(31:60, 2), 36, 'b', 'filled');
title(['Array']);