function A = poolgenU_ffd(meanHD, deltaD1, deltaK1, deltaHD, deltaD2, deltaK2)
% Generates all combinations of sensor parameters
% Inputs:
%   meanHD, deltaD1, deltaK1, deltaHD, deltaD2, deltaK2 - vectors of parameter values
% Output:
%   A - matrix with all K1-K6 combinations and input parameters

    % Input validation
    validateattributes(meanHD, {'numeric'}, {'vector', 'real', 'finite'});
    validateattributes(deltaD1, {'numeric'}, {'vector', 'real', 'finite'});
    validateattributes(deltaK1, {'numeric'}, {'vector', 'real', 'finite'});
    validateattributes(deltaHD, {'numeric'}, {'vector', 'real', 'finite'});
    validateattributes(deltaD2, {'numeric'}, {'vector', 'real', 'finite'});
    validateattributes(deltaK2, {'numeric'}, {'vector', 'real', 'finite'});
    
    % Ensure column vectors
    meanHD = meanHD(:);
    deltaD1 = deltaD1(:);
    deltaK1 = deltaK1(:);
    deltaHD = deltaHD(:);
    deltaD2 = deltaD2(:);
    deltaK2 = deltaK2(:);
    
    % Generate all combinations using ndgrid
    [A_grid, B_grid, C_grid, D_grid, E_grid, F_grid] = ndgrid(meanHD, deltaD1, deltaK1, deltaHD, deltaD2, deltaK2);
    
    % Flatten to column vectors for easier manipulation
    A_vec = A_grid(:);
    B_vec = B_grid(:);
    C_vec = C_grid(:);
    D_vec = D_grid(:);
    E_vec = E_grid(:);
    F_vec = F_grid(:);
    
    % Calculate K values: K1-K6 based on formulas
    % K1 = meanHD - deltaD1/2
    % K2 = K1 + deltaK1 = meanHD - deltaD1/2 + deltaK1
    % K3 = K2 + deltaK1 = meanHD - deltaD1/2 + 2*deltaK1
    % K4 = meanHD + deltaHD/2
    % K5 = K4 + deltaD2 = meanHD + deltaHD/2 + deltaD2
    % K6 = K5 + deltaK2 = meanHD + deltaHD/2 + deltaD2 + deltaK2
    
    K1 = A_vec - D_vec/2;
    K4 = A_vec + D_vec/2;
    K2 = K1 + B_vec;
    K5 = K4 + E_vec;
    K3 = K2 + C_vec;
    K6 = K5 + F_vec;
    
    % Combine K values and input parameters
    M = [A_vec, B_vec, C_vec, D_vec, E_vec, F_vec];
    
    % Create output matrix with K values and parameters
    B_out = [K1, K2, K3, K4, K5, K6, M];
    C_out = [K1, K2, K3, K4, K6, K5, M];
    
    % Concatenate both variants
    A = [B_out; C_out];
    
    % Display summary
    fprintf('Generated %d total combinations\n', size(A,1));
    fprintf('Number of input combinations: %d\n', length(A_vec));
    fprintf('Output matrix size: %d x %d\n', size(A,1), size(A,2));
end
