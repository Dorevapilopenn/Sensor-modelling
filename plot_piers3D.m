% filepath: c:\Users\Propanone\Desktop\Programming\Sensor\Modelling\Sensor-modelling\plot_piers.m
vals = GUI_piers();  % Get values from GUI
C_0 = [vals(3) vals(4) vals(5) vals(6) vals(7)];  % conc A, B1, B2, C respectively
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
title("SC4PUT");

nexttile;
plot(IDA_c{2}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("SC4TYR");

nexttile;
plot(IDA_c{3}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("SC6PUT");

nexttile;
plot(IDA_c{4}, '.-');
legend(lgnd);
xlabel('# Samples'); ylabel('[species]');
title("SC6TYR");

fig1 = figure(2);
% Create a tiled layout (2x3) with compact spacing.
t = tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

nexttile;
plot(lam,A{1},'.-');
legend(lgnd);
xlabel('lambda'); ylabel('ABS');
title("SC4");

D = [IDA_D{1}; IDA_D{2}];
S1trB1 = D(1:50, :);
S1trB2 = D(91:140, :);
S1tB1 = D(51:90, :);
S1tB2 = D(141:180, :);
S1tr  =  [S1trB1, 1*ones(50, 1);
          S1trB2, 2*ones(50, 1)];
S1t   =  [S1tB1, 1*ones(40, 1);
          S1tB2, 2*ones(40, 1)];
S1tr = S1tr(randperm(size(S1tr, 1)), :);  % Shuffle the training set
S1t = S1t(randperm(size(S1t, 1)), :);  % Shuffle the test set
S1trD = S1tr(:, 1:end-1);  % Training data
S1trL = S1tr(:, end);  % Training labels 
S1tD = S1t(:, 1:end-1);  % Test data
S1tL = S1t(:, end);  % Test labels
nexttile;
h1 = plot(lam, D(1:90, :), 'y'); hold on; 
h2 = plot(lam, D(91:180, :), 'g');
legend([h1(1), h2(1)], {'PUT', 'TYR'});
xlabel('lambda'); ylabel('ABS');
title("SC4");

nexttile;
[~,Sc, ~, ~, explained] = pca(D);
plot(Sc(1:90,1),Sc(1:90,2),'r*');hold on; plot(Sc(91:180,1),Sc(91:180,2),'b*')
xlabel("PC1 " + explained(1)); ylabel("PC2 " + explained(2));
title(['Single(SC4)']);

nexttile;
plot(lam,A{2},'.-');
legend(lgnd);
xlabel('lambda'); ylabel('ABS');
title("SC6");

D = [IDA_D{3}; IDA_D{4}];
S2trB1 = D(1:50, :);
S2trB2 = D(91:140, :);
S2tB1 = D(51:90, :);
S2tB2 = D(141:180, :);
S2tr  =  [S2trB1, 1*ones(50, 1);
          S2trB2, 2*ones(50, 1)];
S2t   =  [S2tB1, 1*ones(40, 1);
          S2tB2, 2*ones(40, 1)];
S2tr = S2tr(randperm(size(S2tr, 1)), :);  % Shuffle the training set
S2t = S2t(randperm(size(S2t, 1)), :);  % Shuffle the test set
S2trD = S2tr(:, 1:end-1);  % Training data
S2trL = S2tr(:, end);  % Training labels
S2tD = S2t(:, 1:end-1);  % Test data
S2tL = S2t(:, end);  % Test labels
nexttile;
h1 = plot(lam, D(1:90, :), 'y'); hold on; 
h2 = plot(lam, D(91:180, :), 'g');
legend([h1(1), h2(1)], {'PUT', 'TYR'});
xlabel('lambda'); ylabel('ABS');
title("SC6");

nexttile;
[~,Sc, ~, ~, explained] = pca(D);
plot(Sc(1:90,1),Sc(1:90,2),'r*');hold on; plot(Sc(91:180,1),Sc(91:180,2),'b*')
xlabel("PC1 " + explained(1)); ylabel("PC2 " + explained(2));
title(['Single(SC6)']);

fig2 = figure(3);
D = [IDA_D{1}, IDA_D{3}; IDA_D{2},IDA_D{4}];
AtrB1 = D(1:50, :);
AtrB2 = D(91:140, :);
AtB1 = D(51:90, :);
AtB2 = D(141:180, :);
Atr  =  [AtrB1, 1*ones(50, 1);
          AtrB2, 2*ones(50, 1)];
At   =  [AtB1, 1*ones(40, 1);
          AtB2, 2*ones(40, 1)];
Atr = Atr(randperm(size(Atr, 1)), :);  % Shuffle the training set
At = At(randperm(size(At, 1)), :);  % Shuffle the test set
AtrD = Atr(:, 1:end-1);  % Training data
AtrL = Atr(:, end);  % Training labels
AtD = At(:, 1:end-1);  % Test data
AtL = At(:, end);  % Test labels
[~,Sc, ~, ~, explained] = pca(D);
plot(Sc(1:90,1),Sc(1:90,2),'r*');hold on; plot(Sc(91:180,1),Sc(91:180,2),'b*')
xlabel("PC1 " + explained(1)); ylabel("PC2 " + explained(2));
title(['Array']);