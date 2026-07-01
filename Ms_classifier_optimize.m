function result = Ms_classifier_optimize(varargin)
%Ms_classifier_optimize Sobol/doubling optimizer for single mixed-host sensor.
%
% Default experiment:
%   Sobol variables: HD, deltaD, deltaK, deltaC, deltaG, PN, logCs, r
%   Smoothing grid per nLV: LDA, logistic, and SVM_C = 10.^(-4:0.25:2)
%   Checkpoints: 256, 512, ..., 8192
%   Output JSONs: Ms_optimal_top5.json and Ms_optimal_top1.json
%
% Example smoke run:
%   Ms_classifier_optimize('maxSobolN', 16, 'initialSobolN', 4, ...
%       'stabilityCount', 2, 'nLVGrid', 2:3, 'svmLogCGrid', -1:1);

cfg = parse_options(varargin{:});

checkpoints = cfg.initialSobolN;
while checkpoints(end) < cfg.maxSobolN
    checkpoints(end + 1) = checkpoints(end) * 2; %#ok<AGROW>
end
if checkpoints(end) > cfg.maxSobolN
    checkpoints(end) = cfg.maxSobolN;
end

[nLVByPerm, classifierByPerm, svmLogCByPerm, svmCByPerm] = smoothing_grid(cfg.nLVGrid, cfg.svmLogCGrid);
nPerm = numel(nLVByPerm);
nSamples = cfg.maxSobolN;
[knobs, constants, sampleInfo] = make_idas_sobol_inputs(nSamples, cfg);

if exist(cfg.checkpointFile, 'file')
    S = load(cfg.checkpointFile, 'eff', 'done', ...
        'nLVByPerm', 'classifierByPerm', 'svmLogCByPerm', 'svmCByPerm', 'cfgSaved');
    checkpointCompatible = isfield(S, 'cfgSaved') && ...
        isfield(S.cfgSaved, 'algorithmVersion') && ...
        strcmp(S.cfgSaved.algorithmVersion, cfg.algorithmVersion) && ...
        isfield(S, 'nLVByPerm') && isequal(S.nLVByPerm(:), nLVByPerm(:)) && ...
        isfield(S, 'classifierByPerm') && isequal(S.classifierByPerm(:), classifierByPerm(:)) && ...
        isfield(S, 'svmLogCByPerm') && isequaln(S.svmLogCByPerm(:), svmLogCByPerm(:)) && ...
        isfield(S, 'eff') && size(S.eff, 2) == nPerm && ...
        isfield(S, 'done') && numel(S.done) == size(S.eff, 1);

    if checkpointCompatible
        eff = NaN(nSamples, nPerm, 'single');
        done = false(nSamples, 1);
        reuseRows = min(size(S.eff, 1), nSamples);
        eff(1:reuseRows, :) = S.eff(1:reuseRows, :);
        done(1:reuseRows) = S.done(1:reuseRows);
        cfgSaved = cfg; %#ok<NASGU>
    else
        warning('Ignoring incompatible checkpoint file %s; starting fresh.', cfg.checkpointFile);
        eff = NaN(nSamples, nPerm, 'single');
        done = false(nSamples, 1);
        cfgSaved = cfg; %#ok<NASGU>
    end
else
    eff = NaN(nSamples, nPerm, 'single');
    done = false(nSamples, 1);
    cfgSaved = cfg; %#ok<NASGU>
end

AtrL = [zeros(75, 1); ones(75, 1)];
AtL = [zeros(25, 1); ones(25, 1)];

if cfg.useParallel
    ensure_parallel_pool(cfg.workers);
end

CsVec = sampleInfo.Cs;
C0Vec = sampleInfo.C0;
stopAt = checkpoints(end);

for c = 1:numel(checkpoints)
    checkpointN = checkpoints(c);
    pendingIdx = find(~done(1:checkpointN));

    if ~isempty(pendingIdx)
        fprintf('Evaluating checkpoint chunk to N=%d: %d pending Sobol rows x %d smoothing settings.\n', ...
            checkpointN, numel(pendingIdx), nPerm);

        checkpointTimer = tic;
        overallDoneAtStart = nnz(done);
        if cfg.useParallel
            futures = parallel.FevalFuture.empty(numel(pendingIdx), 0);
            for q = 1:numel(pendingIdx)
                rowIdx = pendingIdx(q);
                futures(q) = parfeval(@evaluate_idas_sobol_row, 3, rowIdx, constants(rowIdx, :), ...
                    CsVec(rowIdx), C0Vec(rowIdx), nLVByPerm, classifierByPerm, svmCByPerm, AtrL, AtL, cfg.baseSeed);
            end

            for q = 1:numel(pendingIdx)
                rowTimer = tic;
                [completedFutureIdx, rowIdx, effRow, workerSeconds] = fetchNext(futures);
                eff(rowIdx, :) = effRow;
                done(rowIdx) = true;
                didSave = false;
                if mod(q, cfg.saveEvery) == 0 || q == numel(pendingIdx)
                    save_checkpoint(cfg.checkpointFile, eff, done, ...
                        nLVByPerm, classifierByPerm, svmLogCByPerm, svmCByPerm, cfgSaved);
                    didSave = true;
                end
                if should_print_progress(q, numel(pendingIdx), didSave, cfg.progressEvery)
                    print_progress(rowIdx, checkpointN, q, numel(pendingIdx), overallDoneAtStart, ...
                        nSamples, checkpointTimer, toc(rowTimer), workerSeconds, nPerm, didSave, cfg.saveEvery);
                end
                futures(completedFutureIdx) = [];
            end
        else
            for q = 1:numel(pendingIdx)
                rowTimer = tic;
                rowIdx = pendingIdx(q);
                [~, effRow, workerSeconds] = evaluate_idas_sobol_row(rowIdx, constants(rowIdx, :), ...
                    CsVec(rowIdx), C0Vec(rowIdx), nLVByPerm, classifierByPerm, svmCByPerm, AtrL, AtL, cfg.baseSeed);
                eff(rowIdx, :) = effRow;
                done(rowIdx) = true;
                didSave = false;
                if mod(q, cfg.saveEvery) == 0 || q == numel(pendingIdx)
                    save_checkpoint(cfg.checkpointFile, eff, done, ...
                        nLVByPerm, classifierByPerm, svmLogCByPerm, svmCByPerm, cfgSaved);
                    didSave = true;
                end
                if should_print_progress(q, numel(pendingIdx), didSave, cfg.progressEvery)
                    print_progress(rowIdx, checkpointN, q, numel(pendingIdx), overallDoneAtStart, ...
                        nSamples, checkpointTimer, toc(rowTimer), workerSeconds, nPerm, didSave, cfg.saveEvery);
                end
            end
        end

        fprintf('Saved checkpoint at Sobol row %d/%d.\n', checkpointN, nSamples);
    else
        fprintf('Checkpoint N=%d already complete in checkpoint file.\n', checkpointN);
    end

    currentCheckpoints = checkpoints(checkpoints <= checkpointN);
    allStable = true;
    for t = 1:numel(cfg.topFractions)
        currentHistory = evaluate_checkpoints(eff, currentCheckpoints, nLVByPerm, ...
            classifierByPerm, svmLogCByPerm, svmCByPerm, cfg.topFractions(t), cfg.stabilityCount);
        allStable = allStable && currentHistory(end).stable;
    end
    if allStable
        stopAt = checkpointN;
        fprintf('Stability reached for all top fractions at Sobol row %d.\n', checkpointN);
        break;
    end
end

completedN = find(done, 1, 'last');
completedCheckpoints = checkpoints(checkpoints <= min(stopAt, completedN));
for t = 1:numel(cfg.topFractions)
    history = evaluate_checkpoints(eff, completedCheckpoints, nLVByPerm, classifierByPerm, ...
        svmLogCByPerm, svmCByPerm, cfg.topFractions(t), cfg.stabilityCount);
    resultForFraction = finalize_result(history, cfg, cfg.topFractions(t), ...
        nLVByPerm, classifierByPerm, svmLogCByPerm, svmCByPerm);
    resultForFraction.modelStats = model_stats_for_json(eff, resultForFraction.selectedAtN);
    write_json(cfg.outputJsons{t}, resultForFraction);
    fprintf('Ms Sobol result written to %s\n', cfg.outputJsons{t});
    if t == 1
        result = resultForFraction;
    else
        result(t) = resultForFraction; %#ok<AGROW>
    end
end
save_checkpoint(cfg.checkpointFile, eff, done, ...
    nLVByPerm, classifierByPerm, svmLogCByPerm, svmCByPerm, cfgSaved);
end

function cfg = parse_options(varargin)
cfg = struct();
cfg.initialSobolN = 256;
cfg.maxSobolN = 262144;
cfg.stabilityCount = 3;
cfg.tolerance = 0.01;
cfg.topFractions = [0.05, 0.01];  % Edit this line: [0.05, 0.01] = top 5% and 1%
cfg.topFraction = [];             % Optional legacy override for a single top fraction.
cfg.nLVGrid = 2:9;
cfg.svmLogCGrid = -4:0.25:2;
cfg.outputJsons = {'Ms_optimal_top5.json', 'Ms_optimal_top1.json'};
cfg.outputJson = '';              % Optional legacy override for a single output file.
cfg.checkpointFile = 'Ms_classifier_optimize_checkpoint.mat';
cfg.saveEvery = 2048;
cfg.progressEvery = 512;
cfg.baseSeed = 1729;
cfg.sobolSkip = 1024;
cfg.sobolLeap = 0;
cfg.useParallel = true;
cfg.workers = 6;
cfg.algorithmVersion = 'Ms_sobol_fast_pls_reuse_svm_squared_hinge_v1';

if mod(numel(varargin), 2) ~= 0
    error('Options must be supplied as name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(cfg, name)
        error('Unknown option "%s".', name);
    end
    cfg.(name) = varargin{k + 1};
end

if ~isempty(cfg.topFraction)
    cfg.topFractions = parse_top_fractions(cfg.topFraction);
else
    cfg.topFractions = parse_top_fractions(cfg.topFractions);
end
if ~isempty(cfg.outputJson)
    cfg.outputJsons = {char(cfg.outputJson)};
end
if ischar(cfg.outputJsons) || isstring(cfg.outputJsons)
    cfg.outputJsons = cellstr(cfg.outputJsons);
end

validateattributes(cfg.initialSobolN, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(cfg.maxSobolN, {'numeric'}, {'scalar', 'integer', '>=', cfg.initialSobolN});
validateattributes(cfg.stabilityCount, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(cfg.tolerance, {'numeric'}, {'scalar', 'positive'});
validateattributes(cfg.topFractions, {'numeric'}, {'vector', '>', 0, '<=', 1});
validateattributes(cfg.nLVGrid, {'numeric'}, {'vector', 'integer', 'positive'});
validateattributes(cfg.svmLogCGrid, {'numeric'}, {'vector', 'real', 'finite'});
validateattributes(cfg.saveEvery, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(cfg.progressEvery, {'numeric'}, {'scalar', 'integer', 'positive'});
if ~iscellstr(cfg.outputJsons) %#ok<ISCLSTR>
    error('outputJsons must be a cell array of character vectors or a string array.');
end
if numel(cfg.outputJsons) ~= numel(cfg.topFractions)
    error('outputJsons must contain one filename per topFractions value.');
end
end

function topFractions = parse_top_fractions(raw)
if iscell(raw)
    topFractions = cellfun(@parse_top_fraction, raw);
elseif isstring(raw) && numel(raw) > 1
    topFractions = arrayfun(@parse_top_fraction, raw);
elseif isnumeric(raw) && numel(raw) > 1
    topFractions = raw;
else
    topFractions = parse_top_fraction(raw);
end
topFractions = topFractions(:).';
end

function topFraction = parse_top_fraction(raw)
if isnumeric(raw)
    topFraction = raw;
elseif isstring(raw) || ischar(raw)
    txt = strtrim(char(raw));
    if ~isempty(txt) && txt(end) == '%'
        topFraction = str2double(strtrim(txt(1:end-1)));
        topFraction = topFraction / 100;
    else
        topFraction = str2double(txt);
        if topFraction > 1
            topFraction = topFraction / 100;
        end
    end
else
    error('topFraction must be numeric or text such as "5%%".');
end

validateattributes(topFraction, {'numeric'}, {'scalar', '>', 0, '<=', 1});
end

function [nLVByPerm, classifierByPerm, svmLogCByPerm, svmCByPerm] = smoothing_grid(nLVGrid, svmLogCGrid)
nLVGrid = nLVGrid(:);
svmLogCGrid = svmLogCGrid(:);
nSvm = numel(svmLogCGrid);
nPerLV = nSvm + 2;
nPerm = numel(nLVGrid) * nPerLV;

nLVByPerm = zeros(nPerm, 1);
classifierByPerm = cell(nPerm, 1);
svmLogCByPerm = NaN(nPerm, 1);
svmCByPerm = NaN(nPerm, 1);

pos = 1;
for i = 1:numel(nLVGrid)
    nLV = nLVGrid(i);

    nLVByPerm(pos) = nLV;
    classifierByPerm{pos} = 'lda';
    pos = pos + 1;

    nLVByPerm(pos) = nLV;
    classifierByPerm{pos} = 'logistic';
    pos = pos + 1;

    idx = pos:(pos + nSvm - 1);
    nLVByPerm(idx) = nLV;
    classifierByPerm(idx) = {'svm-l'};
    svmLogCByPerm(idx) = svmLogCGrid;
    svmCByPerm(idx) = 10 .^ svmLogCGrid;
    pos = pos + nSvm;
end
end

function [knobs, constants, sampleInfo] = make_idas_sobol_inputs(nSamples, cfg)
dim = 8;
if exist('sobolset', 'file') ~= 2
    error('sobolset is required for this experiment.');
end

p = sobolset(dim, 'Skip', cfg.sobolSkip, 'Leap', cfg.sobolLeap);
p = scramble(p, 'MatousekAffineOwen');
u = net(p, nSamples);

ranges = [
     2, 10;   % HD / K1
    -6,  6;   % deltaD
     0, 12;   % deltaK
    -6,  6;   % deltaC
    -8,  8;   % deltaG
     0,  1;   % PN, <0.5=N and >=0.5=P
    -6, -2;   % logCs
    -2,  2];  % r, where C0/Cs = 10^r

scaled = ranges(:, 1).' + u .* (ranges(:, 2).' - ranges(:, 1).');
PN = double(scaled(:, 6) >= 0.5);
knobs = array2table([scaled(:, 1:5), PN, scaled(:, 7:8)], ...
    'VariableNames', {'HD', 'deltaD', 'deltaK', 'deltaC', 'deltaG', 'PN', 'logCs', 'r'});

K1 = scaled(:, 1);
K2 = K1 + scaled(:, 2);
K3 = K2 + scaled(:, 3);
isN = PN == 0;
K2(isN) = K1(isN) + scaled(isN, 2) + scaled(isN, 3);
K3(isN) = K1(isN) + scaled(isN, 2);
[K4, K5] = signed_pair(K1, scaled(:, 4), scaled(:, 5));
constants = [K1, K2, K3, K4, K5];

sampleInfo = table();
sampleInfo.Cs = 10 .^ scaled(:, 7);
sampleInfo.C0 = sampleInfo.Cs .* (10 .^ scaled(:, 8));
end

function [KlowName, KhighName] = signed_pair(baseK, deltaOffset, signedDelta)
lowVal = baseK + deltaOffset;
KlowName = lowVal;
KhighName = lowVal + abs(signedDelta);
neg = signedDelta < 0;
KlowName(neg) = lowVal(neg) - signedDelta(neg);
KhighName(neg) = lowVal(neg);
end

function ensure_parallel_pool(workers)
if exist('parpool', 'file') ~= 2 || exist('gcp', 'file') ~= 2
    error('Parallel execution requires MATLAB Parallel Computing Toolbox.');
end

pool = gcp('nocreate');
if isempty(pool)
    if isempty(workers)
        pool = parpool('local');
    else
        pool = parpool('local', workers);
    end
elseif ~isempty(workers) && pool.NumWorkers ~= workers
    delete(pool);
    pool = parpool('local', workers);
end
set_pool_idle_timeout(pool);
end

function set_pool_idle_timeout(pool)
if isprop(pool, 'IdleTimeout')
    try
        pool.IdleTimeout = Inf;
    catch ME
        warning('Could not disable parallel pool IdleTimeout: %s', ME.message);
    end
end
end

function [rowIdx, effRow, workerSeconds] = evaluate_idas_sobol_row(rowIdx, constantsRow, Cs, C0, ...
    nLVByPerm, classifierByPerm, svmCByPerm, AtrL, AtL, baseSeed)
workerTimer = tic;
rng(baseSeed + rowIdx, 'twister');
[~, f, ~, ~] = Ms_f(constantsRow, Cs, C0);
D = [f{1}; f{2}];
D(~isfinite(D)) = 0;

AtrD = [D(1:75, :); D(101:175, :)];
AtD = [D(76:100, :); D(176:200, :)];

nPerm = numel(nLVByPerm);
effRow = NaN(1, nPerm, 'single');
classes = unique(AtrL);
uniqueNLV = unique(nLVByPerm(:)).';

for nLV = uniqueNLV
    idxForLV = find(nLVByPerm == nLV);
    [scoreTrain, scoreTest, yidxTrain] = fast_plsda_scores(AtrD, AtrL, AtD, nLV, true);

    ldaIdx = idxForLV(strcmp(classifierByPerm(idxForLV), 'lda'));
    if ~isempty(ldaIdx)
        ypred = fast_lda_predict_score(scoreTrain, AtrL, yidxTrain, classes, scoreTest);
        effRow(ldaIdx) = single(efficiency_total(ypred, AtL, classes));
    end

    logisticIdx = idxForLV(strcmp(classifierByPerm(idxForLV), 'logistic'));
    if ~isempty(logisticIdx)
        ypred = fast_logistic_predict_score(scoreTrain, yidxTrain, classes, scoreTest);
        effRow(logisticIdx) = single(efficiency_total(ypred, AtL, classes));
    end

    svmIdx = idxForLV(strcmp(classifierByPerm(idxForLV), 'svm-l'));
    for k = 1:numel(svmIdx)
        p = svmIdx(k);
        ypred = fast_svm_predict_score(scoreTrain, AtrL, classes, scoreTest, svmCByPerm(p));
        effRow(p) = single(efficiency_total(ypred, AtL, classes));
    end
end
workerSeconds = toc(workerTimer);
end

function [scoreTrain, scoreTest, yidx] = fast_plsda_scores(X, y, Xnew, nLV, doScale)
[classes, ~, yidx] = unique(y);
K = numel(classes);
N = size(X, 1);
N1 = size(X, 2);

Y = zeros(N, K);
Y(sub2ind(size(Y), (1:N)', yidx)) = 1;

muX = mean(X, 1);
Xc = X - muX;
if doScale
    sX = std(Xc, 0, 1);
    tol = eps(class(X)) * max(abs(Xc(:)));
    if tol == 0
        tol = eps(class(X));
    end
    sX(sX < tol) = 1;
    sX = max(sX, 1e-10);
    Xc = Xc ./ sX;
else
    sX = ones(1, N1);
end

muY = mean(Y, 1);
Yc = Y - muY;

T = zeros(N, nLV);
P = zeros(N1, nLV);
Q = zeros(K, nLV);
W = zeros(N1, nLV);
E = Xc;
F = Yc;
nActual = nLV;

for h = 1:nLV
    S = E' * F;
    [Uw, sw, Vw] = svd(S, 'econ');
    s = diag(sw);
    tol = max(size(S)) * eps(max(s));
    if sum(s > tol) == 0
        nActual = h - 1;
        break;
    end

    w = Uw(:, 1);
    t = E * w;
    tNorm = norm(t);
    if tNorm < tol
        nActual = h - 1;
        break;
    end

    t = t / tNorm;
    p = E' * t;
    q = F' * t;

    T(:, h) = t;
    P(:, h) = p;
    Q(:, h) = q;
    W(:, h) = w;

    E = E - t * p';
    F = F - t * q';

    if (norm(E, 'fro') / max(norm(Xc, 'fro'), eps)) < tol || ...
            (norm(F, 'fro') / max(norm(Yc, 'fro'), eps)) < tol
        nActual = h;
        break;
    end
end

if nActual < 1
    scoreTrain = zeros(N, 1);
    scoreTest = zeros(size(Xnew, 1), 1);
    return;
end

T = T(:, 1:nActual); %#ok<NASGU>
P = P(:, 1:nActual);
Q = Q(:, 1:nActual);
W = W(:, 1:nActual);

PW = P' * W;
[~, Sp, ~] = svd(PW, 'econ');
tol = max(size(PW)) * eps(max(diag(Sp)));
if sum(diag(Sp) > tol) < size(PW, 1)
    Wstar = W * pinv(PW);
else
    Wstar = W / PW;
end

YhatTrain = Xc * Wstar * Q' + muY;
XnewZ = (Xnew - muX) ./ sX;
YhatTest = XnewZ * Wstar * Q' + muY;

muYhat = mean(YhatTrain, 1);
Yhat0 = YhatTrain - muYhat;
[~, ~, Vpca] = svd(Yhat0, 'econ');
coeff = Vpca(:, 1);

scoreTrain = Yhat0 * coeff;
scoreTest = (YhatTest - muYhat) * coeff;
end

function ypred = fast_lda_predict_score(scoreTrain, y, yidx, classes, scoreTest)
K = numel(classes);
classCounts = accumarray(yidx, 1, [K, 1]);
classMeans = accumarray(yidx, scoreTrain, [K, 1]) ./ classCounts;
centered = scoreTrain - classMeans(yidx);
sigma = (centered' * centered) / max(numel(y) - K, 1);
sigma = sigma + eps(class(scoreTrain)) * max(1, abs(sigma));
invSigmaMu = classMeans ./ sigma;
prior = classCounts(:).' / numel(y);
constant = -0.5 * (classMeans(:).' .* invSigmaMu(:).') + log(max(prior, realmin(class(scoreTrain))));
scores = scoreTest * invSigmaMu(:).' + constant;
[~, idx] = max(scores, [], 2);
ypred = classes(idx);
end

function ypred = fast_logistic_predict_score(scoreTrain, yidx, classes, scoreTest)
X = [ones(size(scoreTrain, 1), 1), scoreTrain];
ybin = double(yidx == 2);
beta = zeros(size(X, 2), 1);
ridge = 1e-8;
ridgeMat = diag([0, 1]);

for iter = 1:30
    eta = max(min(X * beta, 35), -35);
    p = 1 ./ (1 + exp(-eta));
    w = max(p .* (1 - p), eps(class(scoreTrain)));
    grad = X' * (ybin - p) - ridge * ridgeMat * beta;
    H = X' * (X .* w) + ridge * ridgeMat;
    step = H \ grad;
    beta = beta + step;
    if norm(step) <= 1e-8 * (1 + norm(beta))
        break;
    end
end

p2 = 1 ./ (1 + exp(-max(min([ones(size(scoreTest, 1), 1), scoreTest] * beta, 35), -35)));
ypred = classes(ones(size(p2)));
ypred(p2 >= 0.5) = classes(2);
end

function ypred = fast_svm_predict_score(scoreTrain, y, classes, scoreTest, C)
y = y(:);
x = scoreTrain(:);
yt = -ones(size(y));
yt(y == classes(2)) = 1;

w = 0;
b = 0;
activePrev = false(size(yt));
for iter = 1:20
    margin = yt .* (w .* x + b);
    active = margin < 1;
    if ~any(active)
        break;
    end

    xa = x(active);
    ya = yt(active);
    nA = numel(xa);
    H = [1 + 2*C*sum(xa.^2), 2*C*sum(xa); ...
         2*C*sum(xa),          2*C*nA + eps];
    rhs = [2*C*sum(ya .* xa); 2*C*sum(ya)];
    theta = H \ rhs;
    w = theta(1);
    b = theta(2);

    if isequal(active, activePrev)
        break;
    end
    activePrev = active;
end

decision = w .* scoreTest(:) + b;
ypred = repmat(classes(1), size(decision));
ypred(decision >= 0) = classes(2);
end

function effTotal = efficiency_total(ypred, ytrue, classes)
effClass = NaN(numel(classes), 1);
for k = 1:numel(classes)
    cls = classes(k);
    posMask = ytrue == cls;
    negMask = ~posMask;
    sensitivity = sum(ypred(posMask) == cls) / max(sum(posMask), 1);
    specificity = sum(ypred(negMask) ~= cls) / max(sum(negMask), 1);
    effClass(k) = sqrt(sensitivity * specificity);
end
effTotal = mean(effClass(~isnan(effClass)));
end

function print_progress(rowIdx, checkpointN, chunkDone, chunkTotal, overallDoneAtStart, ...
    nSamples, checkpointTimer, waitSeconds, workerSeconds, nPerm, didSave, saveEvery)
elapsed = toc(checkpointTimer);
rowsPerMinute = chunkDone / max(elapsed / 60, eps);
remainingChunk = chunkTotal - chunkDone;
etaMinutes = remainingChunk / max(rowsPerMinute, eps);
overallDone = overallDoneAtStart + chunkDone;
if didSave
    saveText = 'saved';
else
    saveText = sprintf('not saved, saveEvery=%d', saveEvery);
end
fprintf(['row %d done | checkpoint %d: %d/%d (%.1f%%) | overall %d/%d (%.1f%%) | ' ...
    '%.2f rows/min | ETA %.1f min | fetch wait %.2fs | worker %.2fs | %d smoothing fits/row | %s\n'], ...
    rowIdx, checkpointN, chunkDone, chunkTotal, 100 * chunkDone / chunkTotal, ...
    overallDone, nSamples, 100 * overallDone / nSamples, rowsPerMinute, etaMinutes, ...
    waitSeconds, workerSeconds, nPerm, saveText);
end

function tf = should_print_progress(chunkDone, chunkTotal, didSave, progressEvery)
tf = didSave || chunkDone == 1 || chunkDone == chunkTotal || mod(chunkDone, progressEvery) == 0;
end

function history = evaluate_checkpoints(eff, checkpoints, nLVByPerm, classifierByPerm, svmLogCByPerm, ...
    svmCByPerm, topFraction, stabilityCount)
history = repmat(struct( ...
    'N', [], ...
    'selectedIndex', [], ...
    'selectedNLV', [], ...
    'selectedClassifier', [], ...
    'selectedSVMLogC', [], ...
    'selectedSVMC', [], ...
    'meanEfficiency', [], ...
    'standardError', [], ...
    'lowerBound', [], ...
    'maxLowerBound', [], ...
    'candidateCount', [], ...
    'topFraction', [], ...
    'stableRunLength', [], ...
    'stable', false), numel(checkpoints), 1);

selected = NaN(numel(checkpoints), 1);
for c = 1:numel(checkpoints)
    N = checkpoints(c);
    E = double(eff(1:N, :));
    if any(isnan(E(:)))
        error('Efficiency matrix has unfinished rows inside checkpoint N=%d.', N);
    end

    mu = mean(E, 1);
    se = std(E, 0, 1) ./ sqrt(N);
    lb = mu - 2 .* se;
    maxLB = max(lb);
    nTop = max(1, ceil(topFraction * numel(lb)));
    [~, lbOrder] = sort(lb(:), 'descend');
    candidateIdx = lbOrder(1:nTop);

    minNLV = min(nLVByPerm(candidateIdx));
    candidateIdx = candidateIdx(nLVByPerm(candidateIdx) == minNLV);
    priority = classifier_priority(classifierByPerm(candidateIdx));
    svmTieC = svmCByPerm(candidateIdx);
    svmTieC(~strcmp(classifierByPerm(candidateIdx), 'svm-l')) = NaN;
    candidateRows = [priority(:), nan_to_inf(svmTieC(:)), candidateIdx(:)];
    [~, order] = sortrows(candidateRows, [1, 2, 3]);
    pick = candidateIdx(order(1));
    selected(c) = pick;

    runLength = 1;
    for j = c - 1:-1:1
        if selected(j) == pick
            runLength = runLength + 1;
        else
            break;
        end
    end

    history(c).N = N;
    history(c).selectedIndex = pick;
    history(c).selectedNLV = nLVByPerm(pick);
    history(c).selectedClassifier = classifierByPerm{pick};
    history(c).selectedSVMLogC = svmLogCByPerm(pick);
    history(c).selectedSVMC = svmCByPerm(pick);
    history(c).meanEfficiency = mu(pick);
    history(c).standardError = se(pick);
    history(c).lowerBound = lb(pick);
    history(c).maxLowerBound = maxLB;
    history(c).candidateCount = numel(candidateIdx);
    history(c).topFraction = topFraction;
    history(c).stableRunLength = runLength;
    history(c).stable = runLength >= stabilityCount;
end
end

function priority = classifier_priority(classifiers)
priority = zeros(numel(classifiers), 1);
for i = 1:numel(classifiers)
    switch classifiers{i}
        case 'logistic'
            priority(i) = 1;
        case 'lda'
            priority(i) = 2;
        case 'svm-l'
            priority(i) = 3;
        otherwise
            priority(i) = 99;
    end
end
end

function x = nan_to_inf(x)
x(isnan(x)) = Inf;
end

function result = finalize_result(history, cfg, topFraction, nLVByPerm, classifierByPerm, svmLogCByPerm, svmCByPerm)
stableAt = find([history.stable], 1, 'first');
if isempty(stableAt)
    finalIdx = numel(history);
    reachedStability = false;
else
    finalIdx = stableAt;
    reachedStability = true;
end

h = history(finalIdx);
result = struct();
result.sensor = 'Ms';
result.reachedStability = reachedStability;
result.stabilityCount = cfg.stabilityCount;
result.tolerance = cfg.tolerance;
result.topFraction = topFraction;
result.initialSobolN = cfg.initialSobolN;
result.maxSobolN = cfg.maxSobolN;
result.selectedAtN = h.N;
result.optimal = struct( ...
    'nLV', h.selectedNLV, ...
    'classifier', h.selectedClassifier, ...
    'SVM_C', h.selectedSVMC, ...
    'SVM_log10_C', h.selectedSVMLogC, ...
    'permutationIndex', h.selectedIndex);
result.selectedStats = rmfield(h, {'selectedNLV', 'selectedClassifier', 'selectedSVMC', 'selectedSVMLogC'});
result.sobolVariables = {'HD', 'deltaD', 'deltaK', 'deltaC', 'deltaG', 'PN', 'logCs', 'r'};
result.sobolRanges = struct( ...
    'HD', [2, 10], ...
    'deltaD', [-6, 6], ...
    'deltaK', [0, 12], ...
    'deltaC', [-6, 6], ...
    'deltaG', [-8, 8], ...
    'PN', [0, 1], ...
    'logCs', [-6, -2], ...
    'r', [-2, 2]);
result.concentrationTransform = struct('Cs', '10^logCs', 'C0_over_Cs', '10^r', 'C0', 'Cs*10^r');
result.smoothingGrid = struct( ...
    'nLV', cfg.nLVGrid, ...
    'classifiersPerNLV', {{'lda', 'logistic', 'svm-l'}}, ...
    'SVM_log10_C', cfg.svmLogCGrid, ...
    'SVM_C', 10 .^ cfg.svmLogCGrid, ...
    'permutationNLV', nLVByPerm, ...
    'permutationClassifier', {classifierByPerm}, ...
    'permutationSVMLog10C', svmLogCByPerm, ...
    'permutationSVMC', svmCByPerm);
result.history = history;
end

function stats = model_stats_for_json(eff, N)
E = double(eff(1:N, :));
stats = struct();
stats.N = N;
stats.meanEfficiency = mean(E, 1);
stats.sdEfficiency = std(E, 0, 1);
stats.standardError = stats.sdEfficiency ./ sqrt(N);
stats.lowerBound = stats.meanEfficiency - 2 .* stats.standardError;
end

function save_checkpoint(checkpointFile, eff, done, ...
    nLVByPerm, classifierByPerm, svmLogCByPerm, svmCByPerm, cfgSaved)
[checkpointDir, checkpointName, checkpointExt] = fileparts(checkpointFile);
if isempty(checkpointDir)
    checkpointDir = pwd;
end
tmpFile = [tempname(checkpointDir), checkpointExt];
save(tmpFile, 'eff', 'done', ...
    'nLVByPerm', 'classifierByPerm', 'svmLogCByPerm', 'svmCByPerm', 'cfgSaved', '-v7');

lastMessage = '';
for attempt = 1:8
    [ok, message] = movefile(tmpFile, checkpointFile, 'f');
    if ok
        return;
    end

    lastMessage = message;
    pause(0.15 * attempt);
end

warning('Atomic checkpoint replace failed: %s. Falling back to direct save.', lastMessage);
save(checkpointFile, 'eff', 'done', ...
    'nLVByPerm', 'classifierByPerm', 'svmLogCByPerm', 'svmCByPerm', 'cfgSaved', '-v7');
safe_delete_file(tmpFile);
end

function safe_delete_file(filename)
if exist(filename, 'file') ~= 2
    return;
end

warningState = warning('off', 'all');
cleanup = onCleanup(@() warning(warningState));
try
    delete(filename);
catch
end
delete(cleanup);
end

function write_json(filename, result)
txt = jsonencode(result, 'PrettyPrint', true);
fid = fopen(filename, 'w');
if fid < 0
    error('Could not open %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', txt);
delete(cleanup);
end
