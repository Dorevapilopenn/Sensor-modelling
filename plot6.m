param= slider_gui_mixed_U();
names= param(:,1);
values= str2double(param(:, 2));
select= [1,2,4,6,7,9,11,12,14,16,17,19,21,22,24];
Iparam= param(select, :);
Inames= Iparam(:,1);
Ivalues= str2double(Iparam(:, 2));
C_0 = [16.5e-6 1.2e-16 1.3e-16 16.5e-6];  % conc H,G,D respectively
[IDA_D IDA_c Ia] = IDA_function(Ivalues, C_0);
[m_D m_c ma] = mixedfunction(values, C_0);
% Plotting 
% Create a string with names and values
headerText = '';
for i = 1:length(names)
    headerText = [headerText, sprintf('%s: %.2f  ', names{i}, values(i,1))];
end

% Split the header text into two lines
midPoint = ceil(length(headerText) / 2); % Find the midpoint
lineBreakIndex = find(headerText(1:midPoint) == ' ', 1, 'last'); % Break at the last space before midpoint
headerText1 = headerText(1:lineBreakIndex); % First line
headerText2 = headerText(lineBreakIndex+1:end); % Second line

% Combine the two lines with a newline character
headerText = sprintf('%s\n%s', headerText1, headerText2);

fig1=figure(1)
% Create a tiled layout (2x3) with compact spacing.
t = tiledlayout(2,3, 'TileSpacing','Compact','Padding','Compact');

% Add a super title that appears on top of the entire layout.
sgtitle(headerText, 'FontSize', 10);

nexttile;
plot(m_c{1},'.-');
legend(["D", "GD", "HD"]);
xlabel('# Samples');ylabel('[species]');
title("S1G1");
nexttile;
plot(m_c{2},'.-');
legend(["D", "GD", "HD"]);
xlabel('# Samples');ylabel('[species]');
title("S1G2");
nexttile;
D = [m_D{1}; m_D{2}];
[U,S,v] = svd(D,0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Single(S1) '])
nexttile;
plot(m_c{3},'.-');
legend(["D", "GD", "HD"]);
xlabel('# Samples');ylabel('[species]');
title("S2G1");
nexttile;
plot(m_c{4},'.-');
legend(["D", "GD", "HD"]);
xlabel('# Samples');ylabel('[species]');
title("S2G2");
nexttile;
D = [m_D{3}; m_D{4}];
[U,S,v] = svd(D,0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Array '])
%'D:/Artin/sensor/slides/figures/Figure1.png''D:/Artin/sensor/slides/figures/Figure1.png'exportgraphics(fig1, 'D:/Artin/sensor/slides/figures/Figure1.png', 'Resolution', 300);


headerText = '';
for i = 1:length(Inames)
    headerText = [headerText, sprintf('%s: %.2f  ', Inames{i}, values(i,1))];
end

% Split the header text into two lines
midPoint = ceil(length(headerText) / 2); % Find the midpoint
lineBreakIndex = find(headerText(1:midPoint) == ' ', 1, 'last'); % Break at the last space before midpoint
headerText1 = headerText(1:lineBreakIndex); % First line
headerText2 = headerText(lineBreakIndex+1:end); % Second line

% Combine the two lines with a newline character
headerText = sprintf('%s\n%s', headerText1, headerText2);

fig2=figure(2)
% Create a tiled layout (2x3) with compact spacing.
t = tiledlayout(2,3, 'TileSpacing','Compact','Padding','Compact');

% Add a super title that appears on top of the entire layout.
sgtitle(headerText, 'FontSize', 10);

nexttile;
plot(IDA_c{1},'.-');
legend(['D', "HD"]);
xlabel('# Samples');ylabel('[species]');
title("S1G1");
nexttile;
plot(IDA_c{2},'.-');
legend(['D', "HD"]);
xlabel('# Samples');ylabel('[species]');
title("S1G2");
nexttile;
D = [IDA_D{1}; IDA_D{2}];
[U,S,v] = svd(D,0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Single(S1) '])
nexttile;
plot(IDA_c{3},'.-');
legend(['D', "HD"]);
xlabel('# Samples');ylabel('[species]');
title("S2G1");
nexttile;
plot(IDA_c{4},'.-');
legend(['D', "HD"]);
xlabel('# Samples');ylabel('[species]');
title("S2G2");
nexttile;
D = [IDA_D{3}; IDA_D{4}];
[U,S,v] = svd(D,0);
Sc = U(:,1:3)*S(1:3,1:3);
plot3(Sc(1:30,1),Sc(1:30,2),Sc(1:30,3),'r*')
hold on
plot3(Sc(31:60,1),Sc(31:60,2),Sc(31:60,3),'b*')
title(['Array '])
%exportgraphics(fig2, 'D:/Artin/sensor/slides/figures/Figure2.png', 'Resolution', 300);