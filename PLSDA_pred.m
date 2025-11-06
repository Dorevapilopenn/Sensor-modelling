function [ypred, eff_total, eff_class, diagnostics] = plsda_predict(model, Xnew, ytrue)
% Predict & return efficiencies for PLS-DA model with improved performance
% [ypred, eff_total, eff_class, diagnostics] = plsda_predict(model, Xnew, ytrue)
%
% Inputs:
%   model   - PLS-DA model from plsda_train
%   Xnew    - New data matrix (n_new x P)
%   ytrue   - (optional) True class labels for efficiency calculation
%
% Outputs:
%   ypred     - Predicted class labels
%   eff_total - Overall accuracy (if ytrue provided)
%   eff_class - Per-class accuracy [K x 1]
%   diagnostics - Struct with prediction details

% Input validation
if ~isstruct(model) || ~isfield(model, 'Wstar')
    error('Invalid model structure');
end
if ~isnumeric(Xnew) || size(Xnew, 2) ~= size(model.Wstar, 1)
    error('Xnew must be numeric with %d columns', size(model.Wstar, 1));
end

% Center and scale new X data
try
    Xz = (Xnew - model.muX) ./ model.sX;
catch
    error('Error in X preprocessing. Check dimensions and scaling.');
end

% Project using W* to get PLS scores (faster than sequential projection)
Tnew = Xz * model.Wstar;    % Direct projection to latent space

% Compute Y predictions in one step
Yhat_new = Tnew * model.Q' + model.muY;  % Include mean offset in single operation

% Project to discriminant space efficiently
Tsuper = (Yhat_new - model.Yhat_mean_forPCA) * model.pcaCoeff;

% LDA prediction
[ypred, scores] = predict(model.lda, Tsuper);

% Calculate efficiencies if true labels provided
eff_total = [];
eff_class = [];
if nargin >= 3 && ~isempty(ytrue)
    % Vectorized efficiency calculation
    correct = (ypred == ytrue);
    eff_total = mean(correct);
    
    % Vectorized per-class efficiency
    classes = model.classes;
    K = numel(classes);
    eff_class = NaN(K, 1);
    
    % Use logical indexing for faster class-wise calculations
    for k = 1:K
        mask = (ytrue == classes(k));
        if any(mask)
            eff_class(k) = mean(correct(mask));
        end
    end
end

% Pack diagnostics if requested
if nargout >= 4
    diagnostics = struct(...
        'Yhat_new', Yhat_new,...
        'Tplsscores', Tnew,...
        'Tsuper', Tsuper,...
        'lda_scores', scores);
else
    diagnostics = [];
end
end
