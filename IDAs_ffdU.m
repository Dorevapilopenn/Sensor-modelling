HD= 5:1:11;
deltaD= -6:1.5:6;
deltaK= 0:2:14;

A = IDAs_poolgenU_ffd(HD, deltaD, deltaK);

% Save the results
save('D:\Artin\sensor\Sensor-modelling\IDAs_cons9U.mat', 'A');
