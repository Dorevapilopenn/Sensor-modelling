% ffd_a.m — OPTIMIZED: checkpointed, fast, memory-efficient processing

% Parameters
checkpointFile = 'Ma_res9U_checkpoint.mat';
finalOutFile = 'Ma_res9U_sorted.mat';
workers = 6;          % reduced for stability (adjust if you have >16GB RAM)
chunkSize = 500;      % smaller chunks = less RAM per worker batch

% Initialize parallel pool if not already running
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', workers);
end

% LOAD A ONCE: Use matfile for efficient row access (avoid full load)
fprintf('Checking MAT file format...\n');
if exist('Ma_cons9U.mat', 'file')
    try
        % v7.3 MAT files are HDF5-backed and support efficient partial reads.
        h5info('Ma_cons9U.mat');
        disp('MAT file already in v7.3 format.');
    catch
        % Not v7.3, convert it
        fprintf('Converting Ma_cons9U.mat to v7.3 format (more memory efficient)...\n');
        data = load('Ma_cons9U.mat');
        save('Ma_cons9U.mat', '-struct', 'data', '-v7.3');
        fprintf('✓ Conversion complete.\n\n');
    end
else
    error('Ma_cons9U.mat not found!');
end

mf = matfile('Ma_cons9U.mat', 'Writable', false);
numRows = size(mf, 'A', 1);
if numRows < 1
    error('Invalid size of matrix A. Check your data.');
end

% Load or initialize checkpoint
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
    B = zeros(numRows, 3, 'single');  % use single precision to save 50% memory
    lastIdx = 0;
end

% Main loop: process in chunks
AtrL = [zeros(50, 1); ones(50, 1)];
AtL = [zeros(40, 1); ones(40, 1)];
tic;
while lastIdx < numRows
    idx = lastIdx + 1 : min(numRows, lastIdx + chunkSize);
    nIdx = numel(idx);
    
    % Preallocate chunk
    B_chunk = zeros(nIdx, 3, 'single');
    A_chunk = mf.A(idx, :);

    try
        % PARFOR over an in-memory chunk to avoid one HDF5 read per row.
        parfor (kk = 1:nIdx, workers)
            cons = A_chunk(kk, :);

            % Compute response (user function)
            [~, f, ~, ~] = Ma_f(cons, 0.01, 0.0002);

            % Build D and extract regions (no intermediate temp arrays)
            AtrB1 = f{1}(1:50, :);
            AtrB2 = f{3}(1:50, :);
            AtB1  = f{2}(1:40, :);
            AtB2  = f{4}(1:40, :);

            % Construct compact dataset (eliminate unnecessary arrays)
            AtrD = [AtrB1; AtrB2];

            AtD = [AtB1; AtB2];

            % Train and predict
            model = PLSDA_train(AtrD, AtrL, 2, true);
            [~, eff_total, eff_class] = PLSDA_pred(model, AtD, AtL);

            % Store as single to reduce memory
            B_chunk(kk, :) = single([eff_total, eff_class(1), eff_class(2)]);
        end
    catch ME
        % Save current progress before rethrowing
        try
            tmpFile = [checkpointFile '.tmp'];
            save(tmpFile, 'B', 'lastIdx', '-v7.3');
            movefile(tmpFile, checkpointFile);
            fprintf('Saved checkpoint after error: lastIdx = %d / %d (%.1f min)\n', ...
                lastIdx, numRows, toc/60);
        catch
            warning('Failed to save checkpoint during error handling.');
        end
        rethrow(ME);
    end

    % Copy chunk into main array
    B(idx, :) = B_chunk;
    lastIdx = idx(end);

    % Atomic-ish save: write tmp then move
    tmpFile = [checkpointFile '.tmp'];
    save(tmpFile, 'B', 'lastIdx', '-v7.3');
    movefile(tmpFile, checkpointFile);

    % Progress with timing
    elapsed = toc;
    rate = (lastIdx - 0) / (elapsed / 60);  % rows per minute
    remaining = (numRows - lastIdx) / rate;
    fprintf('\x1b[1;32m▶ Progress: %d/%d (%.1f%%) | %.1f min elapsed | ~%.1f min remaining\x1b[0m\n', ...
        lastIdx, numRows, 100*lastIdx/numRows, elapsed/60, remaining);
end

elapsed = toc;
fprintf('Processing complete in %.1f minutes!\n\n', elapsed/60);

% FINAL COMBINE: Keep A compressed to save memory
% Instead of loading full A, use matfile directly for final assembly
clear A_chunk;
fprintf('Assembling final results...\n');

% Write C_sorted directly to file without full in-memory assembly
[~, idx_sorted] = sort(B(:, 1), 'descend');  % sort by efficiency score

% Stream-write to final file to avoid holding full matrix
finalFile_tmp = [finalOutFile '.tmp'];
matfile_out = matfile(finalFile_tmp, 'Writable', true);

[nrows, ncols_A] = size(mf, 'A');
matfile_out.C_sorted(nrows, ncols_A + 3) = single(0);

% Sorted indices are irregular, which matfile cannot use directly.
% Load A once, then stream sorted output chunks to disk.
A_all = mf.A;
for batch_start = 1:5000:nrows
    batch_end = min(batch_start + 4999, nrows);
    batch_idx = idx_sorted(batch_start:batch_end);
    matfile_out.C_sorted(batch_start:batch_end, :) = [A_all(batch_idx, :), B(batch_idx, :)];
end
clear A_all;

% Plot first 100 scores
plot(B(idx_sorted(1:min(100, end)), 1));
xlabel('Rank'); ylabel('Efficiency Score'); grid on;

% Replace original file
movefile(finalFile_tmp, finalOutFile);
fprintf('Final results saved to %s\n', finalOutFile);
