function model = plsda_train(X, y, nLV)
    % X: NxP matrix (rows = samples, columns = variables)
    % y: Nx1 vector with class labels (0/1)
    % nLV: number of latent variables

    % basic validation
    if ~all(ismember(unique(y),[0 1]))
        error('Binary PLS-DA requires y in {0,1}');
    end

    % mean-center + scale
    mu = mean(X);
    sigma = std(X);
    Xz = (X - mu) ./ sigma;

    % PLS regression to class indicator
    [XL,~,XS,~,beta,PCTVAR] = plsregress(Xz,y,nLV);

    % predicted scores & threshold
    yhat = [ones(size(Xz,1),1) Xz] * beta;
    thr = mean(yhat); % naive threshold; tune via CV
    ypred = yhat >= thr;

    model = struct;
    model.mu = mu;
    model.sigma = sigma;
    model.beta = beta;
    model.nLV = nLV;
    model.thr = thr;
    model.PCTVAR = PCTVAR;
end
