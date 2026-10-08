function fit = MRCoSEAManifoldFit(Xk, M, cfg)
%MRCoSEAManifoldFit Local-PCA manifold model of one mode (RM-MEDA style).
%   fit = MRCoSEAManifoldFit(Xk, M, cfg) fits a local PCA to the mode
%   members Xk (n x d, box-normalized decision space). The Pareto set of an
%   M-objective problem is (M-1)-dimensional, so the latent dimension is
%   capped at M-1 (further reduced when the sample is small: n >= 2m+2).
%   Fields of fit:
%     valid   - usability guard: residual variance ratio <= residMax and
%               non-degenerate latent ranges (a thick early cloud or a
%               point-like cluster rejects the model, SBX carries on)
%     mu, V   - mean and eigenvectors (columns, variance-descending)
%     m       - latent dimension; Z = (Xk - mu)*V(:,1:m)
%     zlo/zhi - 5%/95% latent percentiles (the sampling box)
%     sigRes  - per-dimension residual std (sqrt of mean residual eigval)
%     resid   - residual variance ratio 1 - cumvar(m)
%     maxGap  - largest relative gap of the 1-D latent layout (m == 1)
%   (Identical semantics to LMICoSEAManifoldFit; v3.0 fits on the mode
%   members directly -- under mode-relative selection the population
%   already holds clean PS samples of local-PF modes, so the v2.0 archive
%   anchoring is unnecessary.)

    [n, d] = size(Xk);
    fit = struct('valid',false,'mu',zeros(1,d),'V',zeros(d),'m',0, ...
        'Z',zeros(0,0),'zlo',[],'zhi',[],'sigRes',0,'resid',1,'maxGap',NaN);

    m = max(1, min(M - 1, d - 1));
    if n < 2*m + 2
        m = floor((n - 2) / 2);
    end
    if m < 1 || n < 4
        return;                                   % not enough members
    end

    mu = mean(Xk, 1);
    C  = cov(Xk);
    [V, E] = eig((C + C') / 2, 'vector');
    [E, ord] = sort(E, 'descend');
    V = V(:, ord);
    E = max(E, 0);                                % guard tiny negatives
    if numel(E) < d
        return;
    end
    cumv  = cumsum(E) / max(sum(E), realmin);
    % Shrink the latent dimension to the variance spectrum (G4 fix): with
    % few members on a high-dimensional PS the top axes already carry most
    % of the variance, and a full (M-1)-dim latent maximin layout is too
    % sparse to sample meaningfully.
    mVar = find(cumv >= cfg.varThresh, 1);
    if ~isempty(mVar)
        m = min(m, mVar);
    end
    if m < 1
        return;
    end
    resid = 1 - cumv(m);

    Z    = (Xk - mu) * V(:, 1:m);
    qlo  = prctile(Z, 5);
    qhi  = prctile(Z, 95);
    zrng = qhi - qlo;
    if any(~isfinite(zrng)) || all(zrng < 1e-6)
        return;                                   % degenerate latent layout
    end
    if m < d
        sigRes = sqrt(mean(E(m+1:end)));
    else
        sigRes = 0;
    end
    maxGap = NaN;
    if m == 1
        zs = sort(Z);
        maxGap = max(diff(zs)) / max(zrng, realmin);
    end
    fit = struct('valid',resid <= cfg.residMax,'mu',mu,'V',V,'m',m, ...
        'Z',Z,'zlo',qlo,'zhi',qhi,'sigRes',sigRes,'resid',resid, ...
        'maxGap',maxGap);
end
