% Define the base vector
P= [1,5,9,13]

A = IDAa_poolgen_ffd(P);

% Save the results
save('D:\Artin\sensor\Sensor-modelling\IDAa_cons4.mat', 'A');