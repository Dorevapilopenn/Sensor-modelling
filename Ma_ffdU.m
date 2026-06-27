% Define the base vector for parameter space exploration
meanHD= 5:1:11;
deltaD1= -6:1.5:6;
deltaK1= 0:2:14;
deltaC1= -6:1.5:6;
deltaG1= 0:2:14;
deltaHD= -6:1.5:6;
deltaD2= -6:1.5:6;
deltaK2= 0:2:14;
deltaC2= -6:1.5:6;
deltaG2= 0:2:14;

% Generate LHS-sampled parameter combinations (5000 samples from ~36M grid)
A = Ma_poolgenU_ffd(meanHD, deltaD1, deltaK1, deltaC1, deltaG1, deltaHD, deltaD2, deltaK2, deltaC2, deltaG2, 500000);

% Save the results
save('Ma_cons9U.mat', 'A', '-v7.3');
