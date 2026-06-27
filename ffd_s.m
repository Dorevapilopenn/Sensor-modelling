% ffd_s.m — OPTIMIZED: memory-efficient, checkpointed processing

% Initialize parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool('local', 4);  % Reduced to 4 for memory stability
end

% Use matfile for efficient row access (NO full load into memory)
mf = matfile('Ms_cons9U.mat', 'Writable', false);
numRows = size(mf, 'A', 1);
if numRows < 1
    error('Invalid size of matrix A. Check your data.');
end

% Checkpoint & recovery setup
checkpointFile = 'Ms_res9U_checkpoint.mat';
if exist(checkpointFile, 'file')
    S = load(checkpointFile, 'B', 'lastIdx');
    B = S.B;
    lastIdx = S.lastIdx;
    if size(B,1) ~= numRows || size(B,2) < 3
        warning('Checkpoint incompatible — reinitializing.');
        B = zeros(numRows, 3, 'single');
        lastIdx = 0;
    end
else
    B = zeros(numRows, 3, 'single');  % Use single precision (50% RAM savings)
    lastIdx = 0;
end

% Processing with checkpoints
chunkSize = 500;
tic;

while lastIdx < numRows
    idx = lastIdx + 1 : min(numRows, lastIdx + chunkSize);
    nIdx = numel(idx);
    B_chunk = zeros(nIdx, 3, 'single');

    try
        parfor (kk = 1:nIdx, 4)
            % Each worker gets independent RNG stream
            stream = RandStream.getGlobalStream();
            reset(stream, kk + lastIdx);
            
            i = idx(kk);
            
            % Get current combination via matfile (read only needed row)
            cons = mf.A(i, :);
            
            % Calculate spectral response
            [~, f, ~, ~] = Ms_f(cons, 0.01, 0.0002);
            
            % Construct D matrix efficiently
            D = [f{1}; f{2}];
            
            % Check for NaN or Inf values and replace
            if any(isnan(D(:))) || any(isinf(D(:)))
                D(isnan(D) | isinf(D)) = 0;
            end
            
            % Extract submatrices directly
            AtrB1 = D(1:50, :);
            AtrB2 = D(91:140, :);
            AtB1 = D(51:90, :);
            AtB2 = D(141:180, :);
            
            % Construct datasets compactly
            AtrD = [AtrB1; AtrB2];
            AtrL = [zeros(50, 1); ones(50, 1)];
            perm = randperm(100);
            AtrD = AtrD(perm, :);
            AtrL = AtrL(perm);

            AtD = [AtB1; AtB2];
            AtL = [zeros(40, 1); ones(40, 1)];
            perm = randperm(80);
            AtD = AtD(perm, :);
            AtL = AtL(perm);
            
            % Train model and get predictions
            model = PLSDA_train(AtrD, AtrL, 2, true);
            [~, eff_total, eff_class, ~] = PLSDA_pred(model, AtD, AtL);
            
            % Store results as single
            B_chunk(kk, :) = single([eff_total, eff_class(1), eff_class(2)]);
        end
    catch ME
        % Save checkpoint on error
        try
            tmpFile = [checkpointFile '.tmp'];
            save(tmpFile, 'B', 'lastIdx', '-v7.3');
            movefile(tmpFile, checkpointFile);
            fprintf('Checkpoint saved after error: lastIdx = %d / %d\n', lastIdx, numRows);
        catch
            warning('Failed to save checkpoint.');
        end
        rethrow(ME);
    end

    % Update main array and checkpoint
    B(idx, :) = B_chunk;
    lastIdx = idx(end);

    % Atomic save
    tmpFile = [checkpointFile '.tmp'];
    save(tmpFile, 'B', 'lastIdx', '-v7.3');
    movefile(tmpFile, checkpointFile);

    % Progress with ETA
    elapsed = toc;
    rate = lastIdx / (elapsed / 60);  % rows per minute
    remaining = (numRows - lastIdx) / rate;
    fprintf('Progress: %d/%d (%.1f%%) | %.1f min elapsed | ~%.1f min remaining\n', ...
        lastIdx, numRows, 100*lastIdx/numRows, elapsed/60, remaining);
end

elapsed = toc;
fprintf('\nProcessing complete in %.1f minutes!\n\n', elapsed/60);

% Stream final combine without loading full A into memory
fprintf('Assembling final results...\n');

% Sort indices by efficiency score
[~, idx_sorted] = sort(B(:, 1), 'descend');

% Stream-assemble to file
finalFile_tmp = 'Ms_res9U_final.tmp';
matfile_out = matfile(finalFile_tmp, 'Writable', true);

[nrows, ncols_A] = size(mf, 'A');
C_sorted = zeros(nrows, ncols_A + 3, 'single');

% Batch load and combine to avoid memory spike
for batch_start = 1:5000:nrows
    batch_end = min(batch_start + 4999, nrows);
    batch_idx = idx_sorted(batch_start:batch_end);
    C_sorted(batch_start:batch_end, :) = [mf.A(batch_idx, :), B(batch_idx, :)];
end

% Write sorted results
matfile_out.C_sorted = C_sorted;

% Visualization
plot(C_sorted(1:min(100, end), end));
xlabel('Rank'); ylabel('Efficiency Score'); grid on;
title('Ms Results — Top Performers');

% Save with backup
if exist('Ms_res9U.mat', 'file')
    movefile('Ms_res9U.mat', 'Ms_res9U_backup.mat');
end
save('Ms_res9U.mat', 'C_sorted');
fprintf('Results saved to Ms_res9U.mat\n');