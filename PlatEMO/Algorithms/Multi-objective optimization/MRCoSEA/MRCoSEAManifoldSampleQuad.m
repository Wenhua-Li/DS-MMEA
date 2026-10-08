function cand = MRCoSEAManifoldSampleQuad(fq, q, progress, cfg, Problem)
%MRCoSEAManifoldSampleQuad Offspring sampling on a quadratic latent map.
%   cand = MRCoSEAManifoldSampleQuad(fq, q, progress, cfg, Problem) mirrors
%   MRCoSEAManifoldSample verbatim (greedy maximin latent candidates
%   against the members' z coordinates inside the 5-95% box, then
%   reconstruction, residual-subspace Gaussian noise with the budget-decayed
%   multiplier, de-normalization and clipping) except that the
%   reconstruction goes through the quadratic map [z, z^2] * C of
%   MRCoSEAManifoldFitQuad. Non-quadratic fits (fq.quad = false) are
%   delegated to the incumbent sampler so the linear path behaves
%   identically.

    if ~fq.quad
        cand = MRCoSEAManifoldSample(fq, q, progress, cfg, Problem);
        return;
    end
    d  = Problem.D;
    lb = Problem.lower;
    ub = Problem.upper;

    %% Greedy maximin latent candidates (1-D, same algorithm as incumbent)
    nC = min(10*q + 20, 1000);
    Zc = fq.zlo + rand(nC, 1) .* (fq.zhi - fq.zlo);
    dmin = min(pdist2(Zc, fq.Z), [], 2);
    Zsel = zeros(q, 1);
    for t = 1 : q
        [~, b] = max(dmin);
        zNew = Zc(b, :);
        Zsel(t, :) = zNew;
        Zc(b, :) = [];
        dmin(b) = [];
        dmin = min(dmin, pdist2(Zc, zNew));
    end

    %% Reconstruction through the quadratic map + residual noise
    ns = cfg.maniNoise0 + (cfg.maniNoiseEnd - cfg.maniNoise0) * progress;
    Xn = repmat(fq.mu, q, 1) + [Zsel, Zsel.^2 - fq.cz2] * fq.C;
    if fq.m < d
        coef = fq.sigRes * ns * randn(q, d - fq.m);
        Xn = Xn + coef * fq.V(:, fq.m+1:end)';
    end
    X = lb + Xn .* (ub - lb);
    cand = min(max(X, lb), ub);
end
