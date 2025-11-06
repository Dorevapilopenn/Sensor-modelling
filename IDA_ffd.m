% Define the base vector
a = [2, 4, 6];

% Use ndgrid to generate all combinations more efficiently
[i1,i2,i3,i4,i5,i6] = ndgrid(1:3);

% Preallocate matrix A for better performance
A = zeros(3^6, 6);

% Convert indices to values from vector a
A(:,1) = a(i1(:));
A(:,2) = a(i2(:));
A(:,3) = a(i3(:));
A(:,4) = a(i4(:));
A(:,5) = a(i5(:));
A(:,6) = a(i6(:));

% Save the results
save('D:\Artin\sensor\Sensor-modelling\cons3.mat', 'A');
