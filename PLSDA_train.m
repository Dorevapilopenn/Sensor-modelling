function model = PLSDA_train(X, y, nLV, varargin)
% Robust PLS-DA training using PLS2 algorithm with improved numerical stability
% model = PLSDA_train(X, y, nLV)
% model = PLSDA_train(X, y, nLV, classifier)
% model = PLSDA_train(X, y, nLV, classifier, C)
% model = PLSDA_train(X, y, nLV, classifier, C, doScale)
% 
% Inputs:
%   X          - N x P data matrix
%   y          - N x 1 class labels (integers or categorical)
%   nLV        - number of PLS components (recommended <= min(rank(X),K-1))
%   classifier - (optional) 'lda', 'svm-l'/'linear-svm', or 'logistic' (default 'lda')
%   C          - (optional) SVM box constraint for linear SVM (default 1; ignored by logistic)
%   doScale    - (optional) whether to scale X (default true)
%
% Backward compatible: PLSDA_train(X, y, nLV, doScale) still works.

[classifier, C, doScale] = parse_classifier_options(varargin{:});

% Input validation and rank check
if ~isnumeric(X) || ~ismatrix(X)
    error('X must be a numeric matrix');
end

% Remove NaN/Inf values to prevent rank calculation errors
bad_rows = any(isnan(X) | isinf(X), 2);
if any(bad_rows)
    warning('Removing %d rows with NaN/Inf values from X', sum(bad_rows));
    X(bad_rows, :) = [];
    y(bad_rows) = [];
end

[N, N1] = size(X);
if N < 2, error('Need at least 2 samples'); end

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
    if tol == 0
        tol = eps(class(X));
    end
    sX(sX < tol) = 1;  % More robust zero-variance check
    % Clip extremely large values to prevent overflow
    sX = max(sX, 1e-10);  % Avoid division by zero
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
mu_Yhat = mean(Yhat, 1);
Yhat0 = Yhat - mu_Yhat;
[~, ~, Vpca] = svd(Yhat0, 'econ');
coeff = Vpca(:, 1:min(nPC, size(Vpca, 2)));
score = Yhat0 * coeff;

% --- Classifier in PCA space ---
switch classifier
    case 'lda'
        classifierModel = fast_lda_train(score, y, classes, yidx);
        lda = classifierModel;
        svm = [];
        logreg = [];
    case 'svm-l'
        classifierModel = fast_linear_svm_train(score, y, classes, C);
        lda = [];
        svm = classifierModel;
        logreg = [];
    case 'logistic'
        classifierModel = fast_logistic_train(score, yidx, classes);
        lda = [];
        svm = [];
        logreg = classifierModel;
    otherwise
        error('Unsupported classifier "%s". Use "lda", "svm-l", or "logistic".', classifier);
end

% --- Pack model ---
model = struct();
model.classes = classes;
model.nLV = nLV;
model.nPC = nPC;
model.classifier = classifier;
model.svmC = C;
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
model.classifierModel = classifierModel;
model.lda = lda;
model.svm = svm;
model.logreg = logreg;
model.doScale = doScale;
end

function [classifier, C, doScale] = parse_classifier_options(varargin)
nargs = numel(varargin);
classifier = 'lda';
C = [];
doScale = true;

if nargs == 0 || isempty(varargin{1})
    classifier = 'lda';
    C = [];
    doScale = true;
elseif (islogical(varargin{1}) || (isnumeric(varargin{1}) && isscalar(varargin{1}))) && nargs == 1
    doScale = logical(varargin{1});
    classifier = 'lda';
    C = [];
elseif isstring(varargin{1}) || ischar(varargin{1})
    classifier = normalize_classifier_name(varargin{1});
    if nargs < 2 || isempty(varargin{2})
        C = [];
    else
        C = varargin{2};
    end
    if nargs < 3 || isempty(varargin{3})
        doScale = true;
    else
        doScale = varargin{3};
    end
else
    error('classifier must be "lda", "svm-l", or a logical doScale value for old calls.');
end

if strcmp(classifier, 'svm-l')
    if isempty(C)
        C = 1;
    end
    if ~isnumeric(C) || ~isscalar(C) || ~isfinite(C) || C <= 0
        error('Linear SVM parameter C must be a positive finite scalar.');
    end
else
    if ~isempty(C) && (islogical(C) || (isnumeric(C) && isscalar(C))) && nargs == 2
        doScale = logical(C);
    end
    C = [];
end
end

function classifier = normalize_classifier_name(classifier)
classifier = lower(strtrim(char(classifier)));
switch classifier
    case {'lda', 'linear-discriminant', 'linear discriminant'}
        classifier = 'lda';
    case {'svm-l', 'svml', 'linear-svm', 'linear svm', 'svm_linear', 'svmlinear'}
        classifier = 'svm-l';
    case {'logistic', 'logreg', 'logistic-regression', 'logistic regression', 'multinomial-logistic', 'multinomial logistic'}
        classifier = 'logistic';
    otherwise
        error('Unsupported classifier "%s". Use "lda", "svm-l", or "logistic".', classifier);
end
end

function svm = fast_linear_svm_train(score, y, classes, C)
if exist('templateSVM', 'file') ~= 2 || exist('fitcecoc', 'file') ~= 2
    error('Linear SVM classifier requires MATLAB Statistics and Machine Learning Toolbox functions templateSVM and fitcecoc.');
end

template = templateSVM(...
    'KernelFunction', 'linear', ...
    'BoxConstraint', C, ...
    'Standardize', false);
svm = fitcecoc(score, y, ...
    'Learners', template, ...
    'ClassNames', classes, ...
    'Coding', 'onevsone');
end

function logreg = fast_logistic_train(score, yidx, classes)
K = numel(classes);
if K ~= 2
    error('Fast logistic regression currently supports binary classification only.');
end

X = [ones(size(score, 1), 1), score];
ybin = double(yidx == 2);
beta = zeros(size(X, 2), 1);
ridge = 1e-8;
ridgeMat = eye(size(X, 2));
ridgeMat(1, 1) = 0;

for iter = 1:50
    eta = X * beta;
    eta = max(min(eta, 35), -35);
    p = 1 ./ (1 + exp(-eta));
    w = max(p .* (1 - p), eps(class(score)));
    grad = X' * (ybin - p) - ridge * ridgeMat * beta;
    H = X' * (X .* w) + ridge * ridgeMat;

    if rcond(H) < 1e-12
        step = pinv(H) * grad;
    else
        step = H \ grad;
    end

    beta = beta + step;
    if norm(step) <= 1e-8 * (1 + norm(beta))
        break;
    end
end

logreg = struct();
logreg.classes = classes;
logreg.beta = beta;
logreg.model = 'binary-irls';
logreg.iterations = iter;
logreg.ResponseName = 'y';
logreg.ClassNames = classes;
end

function lda = fast_lda_train(score, y, classes, yidx)
K = numel(classes);
[N, nPC] = size(score);
classCounts = accumarray(yidx, 1, [K, 1]);

classMeans = zeros(K, nPC);
for j = 1:nPC
    classMeans(:, j) = accumarray(yidx, score(:, j), [K, 1]) ./ classCounts;
end

centered = score - classMeans(yidx, :);
pooled = centered' * centered;

denom = max(N - K, 1);
sigma = pooled / denom;
tol = eps(class(score)) * max(1, norm(sigma, 'fro'));
sigma = sigma + eye(nPC, class(score)) * tol;

lda = struct();
lda.classes = classes;
lda.prior = classCounts(:).' / N;
lda.mu = classMeans;
lda.sigma = sigma;
lda.invSigmaMuT = sigma \ classMeans.';
lda.constant = -0.5 * sum(classMeans .* lda.invSigmaMuT.', 2).' + log(max(lda.prior, realmin(class(score))));
lda.ResponseName = inputname(2);
lda.ClassNames = classes;
end
