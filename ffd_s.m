% Initialize parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool('local', 6); % Replace 12 with your desired number of workers
end

% Load data
vals = load('Ms_cons9U.mat');
A = vals.A;  % Extract matrix A from the structure
% Get the number of rows explicitly
numRows = size(A, 1);
if ~isscalar(numRows)
    error('Invalid size of matrix A. Check your data.');
end

% Preallocate B matrix for storing efficiencies
B = zeros(numRows, 3);
% Use parfor for parallel processing
parfor i = 1:numRows
    % Get current combination
    cons = A(i, :);
    
    % Calculate IDA response
    [~, f, ~, ~] = Ms_f(cons, 0.01, 0.0002);
    
    % Construct D matrix efficiently
    D = [f{1}; f{2}];
    
    % Check for NaN or Inf values and skip or replace
    if any(isnan(D(:))) || any(isinf(D(:)))
        % Replace NaN/Inf with 0 or skip this iteration
        D(isnan(D) | isinf(D)) = 0;
        warning('NaN/Inf detected in iteration %d; replaced with 0', i);
    end
    
    % Extract submatrices more efficiently
    AtrB1 = D(1:50, :);
    AtrB2 = D(91:140, :);
    AtB1 = D(51:90, :);
    AtB2 = D(141:180, :);
    
    % Construct training and test sets with labels
    Atr = [AtrB1, zeros(50, 1);
        AtrB2, ones(50, 1)];
    At = [AtB1, zeros(40, 1);
        AtB2, ones(40, 1)];
    
    % Random permutation
    Atr = Atr(randperm(100), :);  % 100 is total size of Atr
    At = At(randperm(80), :);     % 80 is total size of At
    
    % Split data and labels
    AtrD = Atr(:, 1:end-1);
    AtrL = Atr(:, end);
    AtD = At(:, 1:end-1);
    AtL = At(:, end);
    
    % Train model and get predictions
    model = PLSDA_train(AtrD, AtrL, 2, true);
    [~, eff_total, eff_class, ~] = PLSDA_pred(model, AtD, AtL);
    
    % Store results
    B(i, :) = [eff_total, eff_class(1), eff_class(2)];
    
    
end

fprintf('\nProcessing complete!\n');

% Combine and sort results
C = [A, B];
C_sorted = sortrows(C, 11, 'descend');
plot(C_sorted(:,11))
% Save results
save('Ms_res9U.mat', 'C');