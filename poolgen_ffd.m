function combinedMatrix = generateSymmetricKfs(P)
% GENERATESYMMETRICKFS Generates matrix of combinations for 6-component sensor
% Format: [H1 | (a,b) | H2 H3 H4] where:
%   H1: single host from P
%   (a,b): unique symmetrical pairs from P
%   H2,H3,H4: single hosts from P
%
% Expected rows: N^5*(N+1)/2 where N = length(P)

    % Input validation and conversion to column vector
    validateattributes(P, {'numeric'}, {'vector', 'real', 'finite'});
    P = P(:);
    N = length(P);
    
    % Pre-calculate array sizes
    numUniquePairs = N * (N + 1) / 2;  % For (a,b) pairs
    numSingleHosts = N^4;               % For H1,H2,H3,H4
    totalRows = numSingleHosts * numUniquePairs;
    
    % Preallocate final matrix
    combinedMatrix = zeros(totalRows, 6, 'like', P);
    
    % Generate unique symmetrical pairs once
    [i, j] = find(triu(ones(N)));
    uniquePairs = [P(i), P(j)];  % All possible (a,b) pairs
    
    % Use block processing for memory efficiency
    blockSize = min(1e6, totalRows);
    numBlocks = ceil(totalRows/blockSize);
    
    % Pre-generate index matrices for single hosts
    [H1_idx, H2_idx, H3_idx, H4_idx] = ndgrid(1:N);
    H_indices = [H1_idx(:), H2_idx(:), H3_idx(:), H4_idx(:)];
    
    for block = 1:numBlocks
        startIdx = (block-1)*blockSize + 1;
        endIdx = min(block*blockSize, totalRows);
        
        % Calculate indices for this block
        [hostIdx, pairIdx] = ind2sub([numSingleHosts, numUniquePairs], (startIdx:endIdx)');
        
        % Set single host columns (1,4,5,6)
        host_block_idx = mod(hostIdx-1, size(H_indices,1)) + 1;
        combinedMatrix(startIdx:endIdx, 1) = P(H_indices(host_block_idx, 1));    % H1
        combinedMatrix(startIdx:endIdx, 4:6) = P(H_indices(host_block_idx, 2:4)); % H2,H3,H4
        
        % Set unique pair columns (2,3)
        combinedMatrix(startIdx:endIdx, 2:3) = uniquePairs(pairIdx, :);
    end
    
    % Verification
    expectedRows = (N^5 * (N+1)) / 2;
    actualRows = size(combinedMatrix, 1);
    
    fprintf('Matrix generated:\n');
    fprintf('  Unique pairs (a,b): %d\n', numUniquePairs);
    fprintf('  Single host combinations: %d\n', numSingleHosts);
    fprintf('  Total rows: %d (expected: %d)\n', actualRows, expectedRows);
    
    if actualRows ~= expectedRows
        warning('MATLAB:RowMismatch', ...
                'Row count mismatch: actual=%d, expected=%d', ...
                actualRows, expectedRows);
    end
end