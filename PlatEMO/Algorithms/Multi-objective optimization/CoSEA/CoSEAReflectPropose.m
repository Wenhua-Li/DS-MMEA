function cand = CoSEAReflectPropose(Population, Div, Problem, nCand, cfg, progress)
%CoSEAReflectPropose Manifold-aligned symmetric reflection proposal (v0.5+).
%   cand = CoSEAReflectPropose(Population, Div, Problem, nCand, cfg, progress)
%   clusters the discovered modes in the box-normalized decision space
%   (connected components with radius epsFrac*sqrt(d), clusters with >=2
%   members only), then reflects sampled members about the box center:
%   x' = lb + ub - x, which lands on centrally symmetric sibling Pareto sets.
%   M4 scheduling: cluster-round-robin quota with a small-cluster minimum
%   (minClusterQuota), source points sampled proportionally to their
%   isolation Div (srcWeight), Gaussian noise std decaying linearly from
%   sigma0 to sigmaEnd over the evaluation budget (noiseDecay).

    X  = Population.decs;
    n  = size(X, 1);
    lb = Problem.lower;
    ub = Problem.upper;
    d  = Problem.D;
    Xn = (X - lb) ./ max(ub - lb, eps);

    %% Cluster the discovered modes (connected components)
    Dm  = pdist2(Xn, Xn);
    thr = cfg.epsFrac * sqrt(d);
    adj = Dm < thr;
    adj(1:n+1:end) = false;
    comp   = conncomp(graph(sparse(adj)));
    nComp  = max(comp);
    validK = [];
    for c = 1 : nComp
        if sum(comp == c) >= 2
            validK(end+1) = c; %#ok<AGROW>
        end
    end

    %% Sample source points: round-robin over clusters (small-cluster guard)
    srcIdx = zeros(1, nCand);
    if isempty(validK)
        srcIdx = randi(n, 1, nCand);
        pos    = nCand;   % FIX 2026-09-27: pos was undefined on this path;
                          % first triggered by scattered high-dim populations
                          % (NMMF5). Behaviour on the normal path unchanged.
    else
        order = validK(randperm(numel(validK)));
        pos   = 0;
        round = 0;
        while pos < nCand
            round = round + 1;
            for c = order
                if pos >= nCand
                    break;
                end
                mem = find(comp == c);
                if round <= numel(mem)
                    pos = pos + 1;
                    if cfg.srcWeight
                        srcIdx(pos) = CoSEARoulette(mem, Div(mem));
                    elseif round <= cfg.minClusterQuota
                        srcIdx(pos) = mem(round);
                    else
                        srcIdx(pos) = mem(randi(numel(mem)));
                    end
                end
            end
            if round > n
                break;                    % safety: all members exhausted
            end
        end
    end
    srcIdx = srcIdx(1:pos);

    %% Reflect about the box center, add decaying noise, clip to the box
    Xr = lb + ub - X(srcIdx, :);
    if cfg.noiseDecay
        sig = cfg.sigma0 + (cfg.sigmaEnd - cfg.sigma0) * progress;
    else
        sig = cfg.sigma0;
    end
    Xr   = Xr + sig .* (ub - lb) .* randn(size(Xr));
    cand = min(max(Xr, lb), ub);
end

function idx = CoSEARoulette(candIdx, w)
% Weighted sampling of one index from candIdx with weights w (Inf-safe).
    if all(~isfinite(w)) || all(w == 0)
        w = ones(size(w));
    else
        mx = max(w(isfinite(w)));
        w(isinf(w)) = 2 * mx;
        w(~isfinite(w)) = 0;
    end
    c = cumsum(w);
    r = rand * c(end);
    idx = candIdx(find(c >= r, 1, 'first'));
end
