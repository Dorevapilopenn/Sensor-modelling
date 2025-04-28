function sepVec = pcaSep(PCA1, PCA2)
    % pcaStrictSeparation computes strict separation metrics between two groups
    % given their PCA score matrices.
    %
    % Inputs:
    %   PCA1 - n1 x m matrix (each row is an observation for group 1)
    %   PCA2 - n2 x m matrix (each row is an observation for group 2)
    %
    % Outputs:
    %   sepVec - 1x5 vector:
    %       sepVec(1): Normalized centroid separation = d / sigma_avg.
    %                 (Strict criterion: ratio >= 5 indicates high separation)
    %       sepVec(2): Mahalanobis distance between group centroids (using pooled covariance).
    %                 (Strict threshold: > 16.3 for 3 dimensions, for example)
    %       sepVec(3): Overlap probability estimated via the Bhattacharyya coefficient.
    %                 (Strict criterion: < 0.001 indicates negligible overlap)
    %       sepVec(4): Wilks' Lambda from MANOVA (lower values indicate better separation,
    %                 e.g. < 0.01)
    %       sepVec(5): Permutation test p-value (using Mahalanobis distance as test statistic;
    %                 p < 0.001 required for strict separation)
    %
    % Note: This function assumes that each input matrix is already in PCA space.
    
    %% 1. Normalized Centroid Separation
    
    n1 = size(PCA1, 1);
    n2 = size(PCA2, 1);
    
    % Compute group centroids (means)
    mu1 = mean(PCA1, 1); % 1 x m vector
    mu2 = mean(PCA2, 1);
    
    % Euclidean distance between centroids
    d = norm(mu1 - mu2);
    
    % Compute average within-group spread (average of standard deviations along each PC)
    std1 = std(PCA1, 0, 1); % 1 x m vector
    std2 = std(PCA2, 0, 1);
    sigma_avg = (mean(std1) + mean(std2)) / 2;
    
    % Normalized centroid separation ratio
    normCentroidSep = d / sigma_avg;
    
    %% 2. Mahalanobis Distance between Centroids
    
    % Compute covariance matrices for each group
    cov1 = cov(PCA1);
    cov2 = cov(PCA2);
    
    % Pooled covariance
    pooledCov = ((n1-1)*cov1 + (n2-1)*cov2) / (n1 + n2 - 2);
    
    % Use pseudo-inverse to be safe if pooledCov is near singular.
    mahalDist = sqrt((mu1 - mu2) * pinv(pooledCov) * (mu1 - mu2)');
    
    %% 3. Overlap Probability via Bhattacharyya Coefficient
    
    % Compute average covariance matrix
    Sigma = (cov1 + cov2) / 2;
    epsilon = 1e-8;  % small regularizer
    detCov1 = det(cov1) + epsilon;
    detCov2 = det(cov2) + epsilon;
    detSigma = det(Sigma) + epsilon;
    
    % Bhattacharyya distance
    dMu = (mu1 - mu2)';  % column vector
    D_B = (1/8) * (mu1 - mu2) * pinv(Sigma) * dMu + (1/2) * log(detSigma / sqrt(detCov1*detCov2));
    
    % Overlap probability estimate (smaller value means less overlap)
    overlapProb = exp(-D_B);
    %% Collect outputs in a vector
    sepVec = [normCentroidSep, mahalDist, overlapProb];
    
end
