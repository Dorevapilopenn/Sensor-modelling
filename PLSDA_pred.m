function ypred = plsda_predict(model, X)
    Xz = (X - model.mu) ./ model.sigma;
    yhat = [ones(size(Xz,1),1) Xz] * model.beta;
    ypred = yhat >= model.thr;
end
