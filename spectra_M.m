function A = piers_spectra()
% Define the wavelength range
lam = 330:1:600;
AD = [.657, 500, 26.27;
    0,    0  , 0;
    .465,500,37;
    .174,360,64.1;
    .481,496,32;
    .12,330,42;
    .48,500,35.5;
    .108,383,34.75;
    .484,497,38;
    .112,351,38;];
A= cell(1,2);
for i = 1:5
    A1 = AD((2*i-1),:);
    A2 = AD((2*i),:);
    S(i,:)= A1(1) * gauss(lam, A1(2), A1(3)) + ...
        A2(1) * gauss(lam, A2(2), A2(3));
end
% Define the spectra for two different sets
r1 = [1 2 4];
r2 = [1 2 5];
r3 = [1 3 4];
r4 = [1 3 5];
A{1} = S(r1,:);
A{2} = S(r2,:);
A{3} = S(r3,:);
A{4} = S(r4,:);
end