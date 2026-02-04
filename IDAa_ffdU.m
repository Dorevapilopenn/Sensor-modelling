% Define the base vector
meanHD= 5:1:11;
deltaD1= -6:1.5:6;
deltaK1= 0:2:14;
deltaHD= -8:2:8;
deltaD2= -6:1.5:6;
deltaK2= 0:2:14;

A = IDAa_poolgenU_ffd(meanHD, deltaD1, deltaK1, deltaHD, deltaD2, deltaK2);

% Save the results
save('D:\Artin\sensor\Sensor-modelling\cons9U.mat', 'A'); 
