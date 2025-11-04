function A = piers_spectra()
    % Define the wavelength range
    lam = 330:1:600;                  
    AD = [.194,382,33.35;
          .485,452,18.28;
          .203,388,37.54;
          .423,456,18.35];
    for i = 1:2
        A1 = AD((2*i-1),:);
        A2 = AD((2*i),:);
        S(i,:)= A1(1) * gauss(lam, A1(2), A1(3)) + ...
              A2(1) * gauss(lam, A2(2), A2(3));
    end
    A= S;
end