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
%   eff_total - Overall efficiency (macro-average of per-class efficiency, where per-class efficiency = sqrt(sensitivity * specificity), if ytrue provided)
%   eff_class - Per-class efficiency (sqrt(sensitivity * specificity)) [K x 1]
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
    classes = model.classes;
    K = numel(classes);
    eff_class = NaN(K, 1);
    
    % Ensure ytrue and ypred are comparable (same type)
    % ...existing code...
    
    for k = 1:K
        cls = classes(k);
        pos_mask = (ytrue == cls);    % true positives + false negatives
        neg_mask = ~pos_mask;         % true negatives + false positives
        
        % Counts for positives
        n_pos = sum(pos_mask);
        if n_pos > 0
            tp = sum(ypred(pos_mask) == cls);
            fn = n_pos - tp;
            sensitivity = tp / (tp + fn);    % TP / (TP + FN)
        else
            sensitivity = NaN;
        end
        
        % Counts for negatives
        n_neg = sum(neg_mask);
        if n_neg > 0
            tn = sum(ypred(neg_mask) ~= cls); % predicted not cls and true not cls
            fp = n_neg - tn;
            specificity = tn / (tn + fp);     % TN / (TN + FP)
        else
            specificity = NaN;
        end
        
        % Per-class efficiency: geometric mean of sensitivity and specificity
        if ~isnan(sensitivity) && ~isnan(specificity)
            eff_class(k) = sqrt(sensitivity * specificity);
        else
            eff_class(k) = NaN;
        end
    end
    
    % Overall efficiency: macro-average of per-class efficiencies (ignore NaNs)
    eff_total = mean(eff_class(~isnan(eff_class)));
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
