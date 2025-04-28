function f = gauss(x,x_max,width)

% Practical Data Analysis in Chemistry
% Marcel Maeder & Yorck-Michael Neuhold
% University of Newcastle, Australia, 2006

% Gaussian of half width <width>, centered around <x_max>

f = exp(-log(2)*((x-x_max)/(width/2)).^2);