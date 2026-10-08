function fq = MRCoSEAManifoldFitQuad(Xk, M, cfg)
%MRCoSEAManifoldFitQuad Quadratic 1-D latent manifold model of one mode.
%   fq = MRCoSEAManifoldFitQuad(Xk, M, cfg) extends the incumbent local PCA
%   with a second-order map, but ONLY on the failure class it targets
%   (latent dimension 1, i.e. curved-curve modes):
%     * m >= 2 or an already clean linear fit (resid <= 1 - varThresh):
%       returns the plain MRCoSEAManifoldFit unchanged (incumbent
%       semantics; samplable by the incumbent sampler);
%     * m == 1 and the linear chord is too coarse (resid > 1 - varThresh):
%       fits a per-dimension quadratic correction on the latent coordinate
%       z by least squares, Xn ~ mu + [z, z^2] * C (closed form,
%       deterministic). Validity bar identical to the segmented arm:
%       residual SS / total SS <= 1 - varThresh, otherwise invalid (SBX
%       carries on, exactly like an incumbent rejection).
%   The returned struct keeps all incumbent fields (mu/V/m/Z/zlo/zhi) so
%   MRCoSEAManifoldSampleQuad can reuse the latent box and the anchor
%   coordinates; extra fields: quad (true/false), C (2 x d map), sigRes
%   (residual std along the orthogonal complement, same noise semantics
%   as the incumbent).
%
%   Motivation (2026-10-01 manifold-upgrade race, dev plan S9): single-bend
%   arcs (e.g. MMF4) need one curved model, not a chain of chords.

    thr = 1 - cfg.varThresh;                   % 0.10 with default config
    fq = MRCoSEAManifoldFit(Xk, M, cfg);       % incumbent fields + gates
    fq.quad = false;
    fq.C = [];
    fq.cz2 = NaN;
    if fq.m ~= 1 || ~isfinite(fq.resid) || fq.resid <= thr
        return;                                % linear path / clean fit
    end
    n  = size(Xk, 1);
    z  = (Xk - fq.mu) * fq.V(:, 1);            % n x 1, mean-centered
    % Basis {z, z^2 - mean(z^2)}: both regressors are mean-zero, so the
    % intercept is absorbed and the no-intercept least squares is
    % consistent (a plain z^2 column has a nonzero mean and would skew the
    % map through the origin -- caught by toy U5d, 2026-10-01).
    cz2 = mean(z.^2);
    Phi = [z, z.^2 - cz2];
    C  = Phi \ (Xk - repmat(fq.mu, n, 1));     % 2 x d, deterministic QR
    E  = Xk - (repmat(fq.mu, n, 1) + Phi * C);
    tot = sum((Xk - repmat(fq.mu, n, 1)).^2, 'all');
    if tot <= 0 || ~isfinite(tot)
        return;                                % degenerate members
    end
    residQ = sum(E.^2, 'all') / tot;
    if residQ > thr
        fq.valid = false;                      % quadratic cannot rescue
        fq.resid = residQ;
        return;
    end
    % Residual std along the orthogonal complement (incumbent noise
    % semantics: sampling noise only in the directions the model ignores)
    d = size(Xk, 2);
    if d > 1
        Eo = E * fq.V(:, 2:d);
        sigQ = sqrt(sum(Eo(:).^2) / (n * (d - 1)));
    else
        sigQ = 0;
    end
    fq.quad  = true;
    fq.valid = true;
    fq.resid = residQ;
    fq.C     = C;
    fq.cz2   = cz2;
    fq.sigRes = sigQ;
end
