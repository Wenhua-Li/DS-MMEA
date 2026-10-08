function lab = MRCoSEAModeCluster(X, Problem, cfg)
%MRCoSEAModeCluster Mode identification in the box-normalized decision space.
%   lab = MRCoSEAModeCluster(X, Problem, cfg) labels every solution with a
%   mode index (1..K, ordered by decreasing size) or 0 (unclustered noise).
%   clusterMode = 'eps': connected components with radius epsFrac*sqrt(d),
%   the same semantics as the clustering inside CoSEAReflectPropose (v0.8);
%   'knn': symmetrized k-nearest-neighbor graph, radius-free, k = knnK or
%   round(sqrt(n)) when empty;
%   'quantile' (v3.1 arm): connected components with a data-adaptive radius
%   thr = radMult * median(1-NN distance). The v3.0 radius is an absolute
%   length that ignores the population state; on high-dimensional problems
%   it chains the whole population into one giant mode (G4/G6) and it lets
%   bridging stragglers fuse developing modes (G2-M4 bistability). The
%   median of the nearest-neighbor distances tracks the population's own
%   scale and is robust to the few bridging/outlier members;
%   'pca' (v3.1 arm): same quantile radius, applied in the whitened
%   principal-component space (cumulative variance varThresh, at least 2
%   axes). Orthogonal PCA preserves pairwise distances, so the gain comes
%   from per-axis variance whitening of anisotropic mode layouts.
%   Clusters with fewer than 2 members are noise. (Identical semantics to
%   LMICoSEAModeCluster on the default branch, new prefix for v3.0.)

    n  = size(X, 1);
    d  = Problem.D;
    lb = Problem.lower;
    ub = Problem.upper;
    Xn = (X - lb) ./ max(ub - lb, eps);

    switch lower(cfg.clusterMode)
        case {'knn', 'knnm'}
            if n < 2
                lab = zeros(1, n);
                return;
            end
            if isempty(cfg.knnK)
                k = max(1, round(sqrt(n)));
            else
                k = max(1, round(cfg.knnK));
            end
            k = min(k, n - 1);
            [~, ord] = sort(pdist2(Xn, Xn), 2);
            nn  = ord(:, 2:k+1);
            ii  = repmat((1:n)', k, 1);
            adj = sparse(ii(:), nn(:), true, n, n);
            if strcmpi(cfg.clusterMode, 'knnm')
                adj = adj & adj.';              % 'knnm' (v3.1 arm): mutual
            else                                % knn, radius-free; the
                adj = adj | adj.';              % intersection kills the
            end                                 % single-linkage bridging
            comp = conncomp(graph(sparse(adj)));
            comp = comp(:)';
        case 'epsq'
            % v3.1 Stage-2 arm: the v0.8 ball-graph chain with a
            % data-driven cap. The absolute radius epsFrac*sqrt(d) chains
            % curve-like modes flawlessly (its looseness PROTECTS sparse
            % developing colonies) but exceeds inter-mode gaps on compact
            % high-D populations and fuses everything into one blob. The
            % cap radMult*Q90(1-NN) only ever SHRINKS the radius: loose
            % (no-op) on low-D curves, ~25x tighter on compact high-D
            % populations. One primitive, zero new fields (radMult and
            % nnQuantile already exist).
            if n < 2
                lab = zeros(1, n);
                return;
            end
            Dm = pdist2(Xn, Xn);
            Dm(1:n+1:end) = inf;
            thr = min(cfg.epsFrac * sqrt(d), ...
                cfg.radMult * prctile(min(Dm, [], 2), cfg.nnQuantile));
            adj = Dm < thr;
            adj(1:n+1:end) = false;
            comp = conncomp(graph(sparse(adj)));
            comp = comp(:)';
        case 'knma'
            % v3.1 Stage-2 arm: mutual-knn cores + coarse absorption.
            % Mutual-knn separates established modes without a radius but
            % shreds sparse developing colonies into noise, whose members
            % then die under the global ranking (Stage-2 first look: rec
            % 0.5-0.67 on IDMPM4T1 / MMF16_l2). Remedy: cores are found at
            % the fine mutual-knn scale and never merge; every small
            % fragment (size < 5% of the population, documented constant)
            % joins its nearest core within the v0.8 absolute radius
            % epsFrac*sqrt(d) -- direct distance only, chains are not
            % followed, so the anti-blob separation survives. Fragments
            % beyond that reach stay noise (discovery candidates).
            if n < 2
                lab = zeros(1, n);
                return;
            end
            if isempty(cfg.knnK)
                k = max(1, round(sqrt(n)));
            else
                k = max(1, round(cfg.knnK));
            end
            k = min(k, n - 1);
            [~, ord] = sort(pdist2(Xn, Xn), 2);
            nn  = ord(:, 2:k+1);
            ii  = repmat((1:n)', k, 1);
            adj = sparse(ii(:), nn(:), true, n, n);
            comp = conncomp(graph(sparse(adj & adj.')));
            comp = comp(:)';
            adoptMin = max(2, round(0.05 * n));
            thrA = cfg.epsFrac * sqrt(d);
            for pass = 1 : 3
                uC    = unique(comp);
                szs   = arrayfun(@(g) sum(comp == g), uC);
                big   = uC(szs >= adoptMin);
                small = uC(szs < adoptMin);
                if isempty(big) || isempty(small)
                    break;
                end
                changed = false;
                for si = 1 : numel(small)
                    idxS  = find(comp == small(si));
                    bestD = inf;
                    bestC = 0;
                    for bi = 1 : numel(big)
                        idxB = find(comp == big(bi));
                        dB = min(pdist2(Xn(idxS,:), Xn(idxB,:)), [], 'all');
                        if dB < bestD
                            bestD = dB;
                            bestC = big(bi);
                        end
                    end
                    if bestD <= thrA
                        comp(idxS) = bestC;
                        changed = true;
                    end
                end
                if ~changed
                    break;
                end
            end
        case 'quantile'
            if n < 2
                lab = zeros(1, n);
                return;
            end
            Dm = pdist2(Xn, Xn);
            Dm(1:n+1:end) = inf;
            thr = cfg.radMult * prctile(min(Dm, [], 2), cfg.nnQuantile);
            adj = Dm < thr;
            adj(1:n+1:end) = false;
            comp = conncomp(graph(sparse(adj)));
            comp = comp(:)';
        case 'pca'
            if n < 2
                lab = zeros(1, n);
                return;
            end
            mu = mean(Xn, 1);
            C  = cov(Xn);
            [V, E] = eig((C + C') / 2, 'vector');
            [E, ord] = sort(E, 'descend');
            V = V(:, ord);
            E = max(E, realmin);
            cumv = cumsum(E) / max(sum(E), realmin);
            m = find(cumv >= cfg.varThresh, 1);
            if isempty(m) || m < 2
                m = min(2, d);
            end
            if any(E(1:m) < 1e-12)
                Xp = Xn;                          % degenerate spectrum:
            else                                  % fall back to raw space
                Xp = (Xn - mu) * (V(:, 1:m) ./ sqrt(E(1:m)).');
            end
            Dm = pdist2(Xp, Xp);
            Dm(1:n+1:end) = inf;
            thr = cfg.radMult * median(min(Dm, [], 2));
            adj = Dm < thr;
            adj(1:n+1:end) = false;
            comp = conncomp(graph(sparse(adj)));
            comp = comp(:)';
        otherwise  % 'eps' (v0.8 semantics)
            Dm  = pdist2(Xn, Xn);
            thr = cfg.epsFrac * sqrt(d);
            adj = Dm < thr;
            adj(1:n+1:end) = false;
            comp = conncomp(graph(sparse(adj)));
            comp = comp(:)';
    end

    lab   = zeros(1, n);
    sizes = accumarray(comp.', 1)';
    uC    = find(sizes >= 2);
    if isempty(uC)
        return;
    end
    [~, ord] = sort(sizes(uC), 'descend');       % deterministic mode order
    uC = uC(ord);
    for k = 1 : numel(uC)
        lab(comp == uC(k)) = k;
    end
end
