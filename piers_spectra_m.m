function A = piers_spectra_m()
    % Define the wavelength range
    lam = 330:1:600;                  
    AD = [.657, 500, 26.27;
          0,    0  , 0;
          .194,382,33.35;
          .485,452,18.28;
          .48,500,35.5;
          .108,383,34.75;
          .465,500,37;
          .174,360,64.1;
          .203,388,37.54;
          .423,456,18.35;
          .230,383,39;
          .400,455,18.35
          .630,506,24.5;
          0,    0  , 0;];
    A= cell(1,2);   
    for i = 1:7
        A1 = AD((2*i-1),:);
        A2 = AD((2*i),:);
        S(i,:)= A1(1) * gauss(lam, A1(2), A1(3)) + ...
              A2(1) * gauss(lam, A2(2), A2(3));
    end
    % Define the spectra for two different sets 
    r1 = [1 2 4 6 7];
    r2 = [1 2 3 5 7];
    A{1} = S(r1,:);
    A{2} = S(r2,:);
end