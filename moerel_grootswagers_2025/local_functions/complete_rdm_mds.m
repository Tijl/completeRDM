function [Y, best_d, stress_vals] = complete_rdm_mds(X, d_range)
    % COMPLETE_RDM_MDS Reconstruct RDM by optimising MDS dimensionality.
    %
    %   [Y, best_d, stress_vals] = COMPLETE_RDM_MDS(X, d_range)
    %   finds the optimal dimensionality in d_range that minimises MDS stress,
    %   imputes missing values, and returns the completed RDM.
    %
    %   Inputs:
    %       X       - Incomplete RDM (NaNs allowed)
    %       d_range - Vector of dimensionalities to search (e.g. 1:20)
    %
    %   Outputs:
    %       Y         - Completed RDM
    %       best_d    - Chosen dimensionality
    %       stress_vals - Stress values for all tested dimensions
    %
    %   Author: Tijl Grootswagers
    %   Date: 22 May 2025
    %   Copyright (c) 2025 Tijl Grootswagers. All rights reserved.

    if nargin < 2
        d_range = 2;
    end

    stress_vals = zeros(size(d_range));
    embeddings = cell(size(d_range));

    % Estimate stress for each dimensionality
    for idx = 1:length(d_range)
        d = d_range(idx);
        try
            [Yd, stress] = mdscale(X, d, 'Start', 'random');
            stress_vals(idx) = stress;
            embeddings{idx} = Yd;
        catch
            stress_vals(idx) = Inf;
            embeddings{idx} = [];
        end
    end

    % Select best dimensionality (lowest stress)
    [~, best_idx] = min(stress_vals);
    best_d = d_range(best_idx);
    embedding = embeddings{best_idx};

    if isempty(embedding)
        fprintf('MDS failed to find a valid embedding for all tested dimensions.\n');
        Y=X;
        return
    end

    % Reconstruct full distance matrix from embedding
    Y = squareform(pdist(embedding));

    % Fill in original known values
    Y(~isnan(X)) = X(~isnan(X));
end
