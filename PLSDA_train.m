function model = plsda_train(X, y, nLV, doScale)
% Robust PLS-DA training using PLS2 algorithm with improved numerical stability
% model = plsda_train(X, y, nLV, doScale)
% 
% Inputs:
%   X       - N x P data matrix
%   y       - N x 1 class labels (integers or categorical)
%   nLV     - number of PLS components (recommended <= min(rank(X),K-1))
%   doScale - (optional) whether to scale X (default true)

if nargin < 4, doScale = true; end

% Input validation and rank check
if ~isnumeric(X) || ~ismatrix(X)
    error('X must be a numeric matrix');
end
[N, N1] = size(X);
if N < 2, error('Need at least 2 samples'); end

% Check rank of X
rankX = rank(X, eps(class(X))*max(size(X))*norm(X,'fro'));
if rankX < min(size(X))
    warning('X is rank deficient (rank=%d). Reducing number of components.', rankX);
    nLV = min(nLV, rankX);
end

% --- Convert labels to dummy matrix Y ---
[classes, ~, yidx] = unique(y);
K = numel(classes);
if K < 2, error('Need at least 2 classes'); end
Y = zeros(N, K);
Y(sub2ind(size(Y), (1:N)', yidx)) = 1;

% --- Center and scale X, center Y ---
muX = mean(X, 1);
Xc = X - muX;
if doScale
    sX = std(Xc, 0, 1);
    tol = eps(class(X))*max(abs(Xc(:)));
    sX(sX < tol) = 1;  % More robust zero-variance check
    Xc = Xc ./ sX;
else
    sX = ones(1, N1);
end

muY = mean(Y, 1);
Yc = Y - muY;

% --- Initialize matrices ---
T = zeros(N, nLV);     % X scores
U = zeros(N, nLV);     % Y scores
P = zeros(N1, nLV);     % X loadings
Q = zeros(K, nLV);     % Y loadings
W = zeros(N1, nLV);     % X weights
B = zeros(nLV, nLV);   % Inner regression coeffs

% --- Modified Kernel PLS2 algorithm with improved stability ---
E = Xc;  % E_0 = X
F = Yc;  % F_0 = Y

for h = 1:nLV
    % Calculate weights using SVD for numerical stability
    S = E' * F;  % Covariance matrix
    [Uw, sw, Vw] = svd(S, 'econ');
    
    % Check singular values
    s = diag(sw);
    tol = max(size(S)) * eps(max(s));
    r = sum(s > tol);
    
    if r == 0
        warning('No significant components found at LV %d', h);
        nLV = h - 1;
        break;
    end
    
    % Extract dominant direction with stability check
    w = Uw(:,1);
    q = Vw(:,1);
    
    % Calculate and check scores
    t = E * w;
    u = F * q;
    
    t_norm = norm(t);
    if t_norm < tol
        warning('Near-zero scores at LV %d - reducing model order', h);
        nLV = h - 1;
        break;
    end
    
    % Normalize scores
    t = t / t_norm;
    
    % Calculate loadings with stability check
    p = E' * t;
    q = F' * t;
    
    % Store results
    T(:,h) = t;
    U(:,h) = u;
    P(:,h) = p;
    Q(:,h) = q;
    W(:,h) = w;
    B(h,h) = u' * t;
    
    % Deflate matrices
    E = E - t * p';
    F = F - t * q';
    
    % Check residuals
    if (norm(E,'fro')/norm(Xc,'fro') < tol) || (norm(F,'fro')/norm(Yc,'fro') < tol)
        warning('Residuals too small at LV %d - stopping', h);
        nLV = h;
        break;
    end
end

% Trim matrices if needed
if nLV < size(T,2)
    T = T(:,1:nLV); U = U(:,1:nLV);
    P = P(:,1:nLV); Q = Q(:,1:nLV);
    W = W(:,1:nLV); B = B(1:nLV,1:nLV);
end

% --- Compute projection matrix W* with improved stability ---
PW = P' * W;
[Up, Sp, Vp] = svd(PW, 'econ');
tol = max(size(PW)) * eps(max(diag(Sp)));
r = sum(diag(Sp) > tol);
if r < size(PW,1)
    warning('P''W is rank deficient. Using pseudoinverse for W*.');
    Wstar = W * pinv(PW);
else
    Wstar = W / PW;
end

% --- Compute training predictions ---
Tfull = Xc * Wstar;
Yhat_c = Tfull * Q';
Yhat = Yhat_c + muY;

% --- PCA on predictions for discriminant space ---
nPC = max(1, K-1);
[coeff, score, ~, ~, ~, mu_Yhat] = pca(Yhat, 'Centered', true, ...
    'NumComponents', nPC);

% --- LDA in PCA space ---
lda = fitcdiscr(score, y, 'DiscrimType', 'linear');

% --- Pack model ---
model = struct();
model.classes = classes;
model.nLV = nLV;
model.nPC = nPC;
model.muX = muX;
model.sX = sX;
model.muY = muY;
model.W = W;
model.P = P;
model.Q = Q;
model.B = B;
model.T = T;
model.U = U;
model.Wstar = Wstar;
model.Yhat_train = Yhat;
model.Yhat_mean_forPCA = mu_Yhat;
model.pcaCoeff = coeff;
model.lda = lda;
model.doScale = doScale;
end
