% Define the base vector
HD= 5:1:11;
deltaD= -6:1.5:6;
deltaK= 0:2:14;
deltaC= -6:1.5:6;
deltaG= 0:2:14;

A = Ms_poolgenU_ffd(HD, deltaD, deltaK, deltaC, deltaG);

% Save the results
save('D:\Artin\sensor\Sensor-modelling\Ms_cons9U.mat', 'A');