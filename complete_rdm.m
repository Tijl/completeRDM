function Y = complete_rdm(X, verbose)
    %COMPLETE_RDM Estimate missing values in a representational dissimilarity matrix (RDM).
    %
    %   Y = COMPLETE_RDM(X) returns a completed matrix Y by estimating the 
    %   missing entries (NaNs) in the input matrix X using known distances 
    %   via the Pythagorean theorem.
    %
    %   Y = COMPLETE_RDM(X, VERBOSE) enables optional verbose output. If VERBOSE 
    %   is true, the function prints progress and estimation details during 
    %   reconstruction.
    %
    %   Assumes X is a symmetric distance matrix with NaNs for missing values. 
    %   Estimates are based on shared known distances to other entries. If no 
    %   references are found in a pass, the function will retry in subsequent passes 
    %   until all missing entries are estimated. Also assumes that there are no entire 
    %   rows with missing values.
    %
    %   Inputs:
    %     X        - Incomplete symmetric distance matrix containing NaNs
    %     VERBOSE  - (Optional) flag to print estimation progress (default: false)
    %
    %   Output:
    %     Y        - Completed symmetric distance matrix with all entries filled
    %
    %   Copyright (c) 2025 Tijl Grootswagers
    %
    if nargin < 2
        verbose = false;
    end
    assert(all(sum(~isnan(X))>1),'Cannot reconstruct when a full row is missing')
    Y = X;
    while 1
        [missing_i, missing_j] = find(isnan(triu(Y,1)));
        nmissing = numel(missing_i);
        if ~nmissing % keep going until all missing entries are filled
            return
        end
        for n = 1:nmissing
            i = missing_i(n);
            j = missing_j(n);
            ai = Y(i,:);
            bj = Y(j,:);
            % find reference entries
            known = ~isnan(ai) & ~isnan(bj) & (1:length(ai) ~= i) & (1:length(ai) ~= j);
            if any(known) % if no references, we'll try again in another pass
                a = ai(known);
                b = bj(known);
                % we estimate the missing distance using Pythagorean theorem, 
                % either as the square root of the sum of squares 
                % or the difference of squares of these known distances
                d1 = sqrt(a.^2 + b.^2);
                d2 = sqrt(abs(a.^2 - b.^2));
                % we take the median of all estimates
                d_est = median([d1 d2], 'all');
                % assign estimate to matrix
                Y(i,j) = d_est;
                Y(j,i) = d_est;
                if verbose
                    fprintf('%i/%i Estimated (%d,%d) and (%d,%d) with %.4f using %d references\n', ...
                            n, nmissing, i, j, j, i, d_est, numel(a));
                end
            else
                if verbose
                    fprintf('%i/%i (%d,%d) and (%d,%d) skipping: no references\n', ...
                            n, nmissing, i, j, j, i);
                end
            end
        end
    end
end
