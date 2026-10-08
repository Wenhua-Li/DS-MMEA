function cand = MRCoSEAManifoldSample(fit, q, progress, cfg, Problem)
%MRCoSEAManifoldSample Direct offspring sampling on the learned manifold.
%   cand = MRCoSEAManifoldSample(fit, q, progress, cfg, Problem) draws q
%   latent points by greedy maximin selection against the members' latent
%   coordinates (fills coverage gaps along the PS manifold, targeting
%   IGDX/CR uniformity), reconstructs decision vectors in the latent
%   subspace, and adds residual-subspace Gaussian noise scaled by the
%   local-PCA residual std times a budget-decayed multiplier. The results
%   are de-normalized and clipped to the box.
%   (Identical semantics to LMICoSEAManifoldSample, new prefix for v3.0.)

    m  = fit.m;
    d  = Problem.D;
    lb = Problem.lower;
    ub = Problem.upper;

    %% Greedy maximin latent candidates (fill the largest gaps first).
    %    Incremental update: only the distance to the newest anchor is
    %    computed per round (O(nC) instead of O(nC*nAnchors)).
    nC = min(10*q + 20, 1000);
    Zc = fit.zlo + rand(nC, m) .* (fit.zhi - fit.zlo);
    dmin = min(pdist2(Zc, fit.Z), [], 2);
    Zsel = zeros(q, m);
    for t = 1 : q
        [~, b] = max(dmin);
        zNew = Zc(b, :);
        Zsel(t, :) = zNew;
        Zc(b, :) = [];
        dmin(b) = [];
        dmin = min(dmin, pdist2(Zc, zNew));
    end

    %% Reconstruction with residual-subspace noise
    ns = cfg.maniNoise0 + (cfg.maniNoiseEnd - cfg.maniNoise0) * progress;
    Xn = repmat(fit.mu, q, 1) + Zsel * fit.V(:, 1:m)';
    if fit.m < d
        coef = fit.sigRes * ns * randn(q, d - fit.m);
        Xn = Xn + coef * fit.V(:, fit.m+1:end)';
    end
    X = lb + Xn .* (ub - lb);
    cand = min(max(X, lb), ub);
end
