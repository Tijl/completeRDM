function Y = complete_rdm_graph(X)
    % Store original NaN mask
    mask = isnan(X);
    
    % Replace NaNs with Inf for computation
    G = X;
    G(mask) = Inf;
    n = size(G, 1);
    
    % Multi-pass update for originally-NaN entries only
    changed = true;
    while changed
        changed = false;
        for i = 1:n
            for j = (i+1):n
                if mask(i,j)
                    old_val = G(i,j);
                    for k = 1:n
                        if k ~= i && k ~= j
                            G(i,j) = min(G(i,j), G(i,k) + G(k,j));
                        end
                    end
                    if G(i,j) < old_val
                        changed = true;
                    end
                end
            end
        end
    end

    % Final clean-up
    G(isinf(G)) = NaN;              % Restore NaNs for unreachable pairs
    G(logical(eye(n))) = 0;         % Zero diagonal
    G = triu(G,1) + triu(G,1)';     % Reflect to lower triangle

    Y = G;
end