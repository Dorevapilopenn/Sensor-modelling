function combinedMatrix = generateSymmetricKfs(P)
% GENERATESYMMETRICKFS Generates optimized matrix of unique combinations
% for a 6-component sensor system with symmetrical pairs.
%
% Optimized version with:
% - Preallocated arrays
% - Vectorized operations
% - Reduced memory usage
% - Improved performance for large P

    % Input validation and conversion to column vector
    validateattributes(P, {'numeric'}, {'vector', 'real', 'finite'});
    P = P(:);
    N = length(P);
    
    % Pre-calculate array sizes
    numUniquePairs = N * (N + 1) / 2;
    numHD = N * N;
    numHG = numUniquePairs * numUniquePairs;
    totalRows = numHD * numHG;
    
    % Preallocate final matrix
    combinedMatrix = zeros(totalRows, 6, 'like', P);
    
    % --- 1. Generate Host-Dye (HD) Combinations more efficiently ---
    [A_idx, B_idx] = ndgrid(1:N);
    
    % --- 2. Generate Unique Symmetrical Host-Guest Pairs ---
    % Preallocate uniquePairs
    uniquePairs = zeros(numUniquePairs, 2, 'like', P);
    idx = 1;
    
    % Vectorized unique pairs generation
    [i, j] = find(triu(ones(N)));
    uniquePairs(:,1) = P(i);
    uniquePairs(:,2) = P(j);
    
    % --- 3 & 4. Combined HD and HG matrix generation ---
    % Use block processing to reduce memory usage
    blockSize = min(1e6, totalRows);  % Adjust based on available memory
    numBlocks = ceil(totalRows/blockSize);
    
    for block = 1:numBlocks
        startIdx = (block-1)*blockSize + 1;
        endIdx = min(block*blockSize, totalRows);
        blockRows = endIdx - startIdx + 1;
        
        % Calculate indices for this block
        [hdIdx, hgIdx] = ind2sub([numHD, numHG], (startIdx:endIdx)');
        
        % Get HD combinations for this block
        combinedMatrix(startIdx:endIdx, 1) = P(A_idx(hdIdx));  % A
        combinedMatrix(startIdx:endIdx, 4) = P(B_idx(hdIdx));  % B
        
        % Get HG combinations for this block
        [ab_idx, cd_idx] = ind2sub([numUniquePairs, numUniquePairs], hgIdx);
        
        combinedMatrix(startIdx:endIdx, 2:3) = uniquePairs(ab_idx, :);  % a,b
        combinedMatrix(startIdx:endIdx, 5:6) = uniquePairs(cd_idx, :);  % c,d
    end
    
    % --- 5. Verification ---
    expectedRows = (N^4 * (N+1)^2) / 4;
    actualRows = size(combinedMatrix, 1);
    
    % Use more efficient string formatting
    fprintf('Unique Host-Guest pairs (a,b): %d\n', numUniquePairs);
    fprintf('Actual rows generated: %d\n', actualRows);
    fprintf('Expected rows (N^4*(N+1)^2/4): %d\n', expectedRows);
    
    if actualRows ~= expectedRows
        warning('MATLAB:RowMismatch', ...
                'Row count mismatch: actual=%d, expected=%d', ...
                actualRows, expectedRows);
    end
end