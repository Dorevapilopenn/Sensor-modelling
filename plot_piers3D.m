% filepath: c:\Users\Propanone\Desktop\Programming\Sensor\Modelling\Sensor-modelling\plot_piers.m
vals = GUI_piers();  % Get values from GUI
C_0 = [vals(3) vals(4) vals(5) vals(6)];  % conc A, B1, B2, C respectively
[C_eq, IDA_D, IDA_c, Ia,  A] = IDA_f_piers(C_0, vals(1), vals(2))
lgnd= ["C", "CH+", "AC", "ACH+"];

lam= 330:1:600;  % Define the wavelength range
% Plotting

fig1 = figure(1);
% Create a tiled layout (2x3) with compact spacing.
t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

nexttile;
plot(IDA_c{1}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("SC4B1");

nexttile;
plot(IDA_c{2}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("SC4B2");

nexttile;
plot(IDA_c{3}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("SC6B1");

nexttile;
plot(IDA_c{4}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("SC6B2");

fig1 = figure(2);
% Create a tiled layout (2x3) with compact spacing.
t = tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

nexttile;
plot(lam,A{1},'.-');
legend(lgnd);
xlabel('lambda'); ylabel('ABS');
title("SC4");

D = [IDA_D{1}; IDA_D{2}];
nexttile;
plot(lam, D, 'y');
xlabel('lambda'); ylabel('ABS');
title("SC4");

nexttile;
[~,Sc, ~, ~, explained] = pca(D);
plot(Sc(1:30,1),Sc(1:30,2),'r*');hold on; plot(Sc(31:60,1),Sc(31:60,2),'b*')
xlabel("PC1 " + explained(1)); ylabel("PC2 " + explained(2));
title(['Single(SC4)']);

nexttile;
plot(lam,A{2},'.-');
legend(lgnd);
xlabel('lambda'); ylabel('ABS');
title("SC6");

D = [IDA_D{3}; IDA_D{4}];
nexttile;
plot(lam, D, 'y');
xlabel('lambda'); ylabel('ABS');
title("SC6");

nexttile;
[~,Sc, ~, ~, explained] = pca(D);
plot(Sc(1:30,1),Sc(1:30,2),'r*');hold on; plot(Sc(31:60,1),Sc(31:60,2),'b*')
xlabel("PC1 " + explained(1)); ylabel("PC2 " + explained(2));
title(['Single(SC6)']);

fig2 = figure(3);
D = [IDA_D{1}, IDA_D{3}; IDA_D{2},IDA_D{4}];
[~,Sc, ~, ~, explained] = pca(D);
plot(Sc(1:30,1),Sc(1:30,2),'r*');hold on; plot(Sc(31:60,1),Sc(31:60,2),'b*')
xlabel("PC1 " + explained(1)); ylabel("PC2 " + explained(2));
title(['Array']);