% Define the base vector
meanHD= 5:1:11;
deltaD1= -6:1.5:6;
deltaK1= 0:2:14;
deltaC1= -6:1.5:6
deltaG1= 0:2:14
deltaHD= -8:2:8;
deltaD2= -6:1.5:6;
deltaK2= 0:2:14;
deltaC1= -6:1.5:6
deltaG1= 0:2:14

A = Ma_poolgenU_ffd(meanHD, deltaD1, deltaK1, deltaC1, deltaG1, deltaHD, deltaD2, deltaK2, deltaC1, deltaG1);

% Save the results
save('D:\Artin\sensor\Sensor-modelling\Ma_cons9U.mat', 'A');