function A = Ma_poolgenU_ffd(meanHD, deltaD1, deltaK1, deltaC1, deltaG1, deltaHD, deltaD2, deltaK2, deltaC2, deltaG2, num_lhs_samples)
    % Ma_poolgenU_ffd: Generate parameter combinations using Latin Hypercube Sampling
    % Inputs:
    %   meanHD, deltaD1, deltaK1, deltaC1, deltaG1, deltaHD, deltaD2, deltaK2, deltaC2, deltaG2 - parameter vectors
    %   num_lhs_samples - number of LHS samples (default: 5000)
    % Output:
    %   A - matrix with all K-value combinations and input parameters via LHS sampling
    
    if nargin < 11 || isempty(num_lhs_samples)
        num_lhs_samples = 5000;
    end
    
    % Ensure column vectors
    meanHD = meanHD(:);
    deltaD1 = deltaD1(:);
    deltaK1 = deltaK1(:);
    deltaC1 = deltaC1(:);
    deltaG1 = deltaG1(:);
    deltaHD = deltaHD(:);
    deltaD2 = deltaD2(:);
    deltaK2 = deltaK2(:);
    deltaC2 = deltaC2(:);
    deltaG2 = deltaG2(:);
    
    % Parameter info
    params = {meanHD, deltaD1, deltaK1, deltaC1, deltaG1, deltaHD, deltaD2, deltaK2, deltaC2, deltaG2};
    param_ranges = cellfun(@(p) [min(p), max(p)], params, 'UniformOutput', false);
    grid_size = prod(cellfun(@length, params));
    
    fprintf('=== Ma_poolgenU_ffd: Latin Hypercube Sampling ===\n');
    fprintf('Full grid size: %.2e points\n', grid_size);
    fprintf('LHS samples: %d\n', num_lhs_samples);
    fprintf('Sampling ratio: %.4f%%\n', 100*num_lhs_samples/grid_size);
    
    % Generate Latin Hypercube Samples in normalized [0,1] space
    rng(42);  % Reproducibility
    n_dims = length(params);
    lhs_norm = lhsdesign(num_lhs_samples, n_dims, 'smooth', 'on');
    
    % Convert LHS samples to actual parameter values using nearest neighbor in parameter vectors
    sampled = zeros(num_lhs_samples, n_dims);
    for d = 1:n_dims
        param_vals = params{d};
        % Map [0,1] to parameter indices
        idx = round(lhs_norm(:, d) * (length(param_vals) - 1)) + 1;
        idx = max(1, min(length(param_vals), idx));
        sampled(:, d) = param_vals(idx);
    end
    clear lhs_norm;  % Free memory
    
    % Extract sampled parameter vectors
    A_vec = sampled(:, 1);   % meanHD
    B_vec = sampled(:, 2);   % deltaD1
    C_vec = sampled(:, 3);   % deltaK1
    D_vec = sampled(:, 4);   % deltaC1
    E_vec = sampled(:, 5);   % deltaG1
    F_vec = sampled(:, 6);   % deltaHD
    G_vec = sampled(:, 7);   % deltaD2
    H_vec = sampled(:, 8);   % deltaK2
    I_vec = sampled(:, 9);   % deltaC2
    J_vec = sampled(:, 10);  % deltaG2
    clear sampled;
    
    % Calculate derived K values based on ffd_a logic
    % K-values represent frequency/stiffness parameters
    K1 = A_vec - F_vec/2;
    K6 = A_vec + F_vec/2;
    K2 = K1 + B_vec;
    K3 = K2 + C_vec;
    K4 = K1 + D_vec;
    K5 = K4 + E_vec;
    K7 = K6 + G_vec;
    K8 = K7 + H_vec;
    K9 = K6 + I_vec;
    K10 = K9 + J_vec;
    
    % Store input parameters in M
    M = [A_vec, B_vec, C_vec, D_vec, E_vec, F_vec, G_vec, H_vec, I_vec, J_vec];
    
    % Generate FFD variants (8 permutations as per original design)
    % Each variant represents different orderings of K-values
    variants = cell(8, 1);
    variants{1} = [K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, M];  % B
    variants{2} = [K1, K2, K3, K4, K5, K6, K7, K8, K10, K9, M]; % C
    variants{3} = [K1, K2, K3, K4, K5, K6, K8, K7, K9, K10, M]; % D
    variants{4} = [K1, K2, K3, K4, K5, K6, K8, K7, K10, K9, M]; % E
    variants{5} = [K1, K2, K3, K5, K4, K6, K7, K8, K9, K10, M]; % F
    variants{6} = [K1, K2, K3, K5, K4, K6, K7, K8, K10, K9, M]; % G
    variants{7} = [K1, K2, K3, K5, K4, K6, K8, K7, K9, K10, M]; % H
    variants{8} = [K1, K2, K3, K5, K4, K6, K8, K7, K10, K9, M]; % I
    
    % Clear K values to save memory
    clear K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 M A_vec B_vec C_vec D_vec E_vec F_vec G_vec H_vec I_vec J_vec;
    
    % Deduplicate each variant, then combine
    fprintf('\nProcessing variants:\n');
    A = [];
    variant_names = {'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'};
    
    for v = 1:8
        % Deduplicate this variant
        [variant_unique, ~, ~] = unique(variants{v}, 'rows', 'stable');
        A = [A; variant_unique];
        fprintf('  Variant %s: %d unique rows → total: %d\n', variant_names{v}, size(variant_unique, 1), size(A, 1));
        clear variants{v} variant_unique;
    end
    
    % Final global deduplication
    fprintf('\nFinal deduplication:\n');
    fprintf('  Before: %d rows\n', size(A, 1));
    [A, ~, ~] = unique(A, 'rows', 'stable');
    fprintf('  After: %d rows × %d columns\n', size(A, 1), size(A, 2));
    fprintf('\n=== Complete ===\n');
end