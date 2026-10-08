function segs = MRCoSEAManifoldFitSeg(Xk, M, minS, cfg)
%MRCoSEAManifoldFitSeg Segmented-chord manifold model of one mode.
%   segs = MRCoSEAManifoldFitSeg(Xk, M, minS, cfg) fits ONE mode (n x d,
%   box-normalized decision space) as a small set of linear pieces instead
%   of the incumbent single local PCA:
%     * a piece is accepted when its plain MRCoSEAManifoldFit explains at
%       least varThresh of the variance (resid <= 1 - varThresh; the dual
%       reading of the existing latent-dimension rule -- the model only
%       trusts what it can explain);
%     * otherwise the members are split at the MEDIAN of the 1st principal
%       coordinate (deterministic, closed form) and both halves are
%       recursed;
%     * a piece that cannot be modeled and is too small to split is
%       discarded (SBX carries that territory).
%   The MODE-level maturity gate (caller minS, minModeSize-based) has
%   already accepted the whole mode; a PIECE only needs enough points for
%   the fit itself, so the recursion floor is the statistical minimum of
%   MRCoSEAManifoldFit (2m+2), not minS. Zero new tunable parameters: the
%   split bar reuses varThresh and the piece cap (16) is an implementation
%   constant like the "keep GA alive" N-2 rule. Linear modes return the
%   identical single-segment fit of the incumbent.
%   Fields of one segment: fit (a plain MRCoSEAManifoldFit result,
%   samplable by MRCoSEAManifoldSample verbatim) and size (members).
%
%   Motivation (2026-10-01 manifold-upgrade race, dev plan S9): the
%   incumbent single chord per mode leaves the channel at zero net gain on
%   the 13 curved-PS problems of the benchmark (truth-mode residual >
%   0.20); piecewise chords target exactly that failure class.

    SEG_CAP = 16;                              % implementation constant
    thr  = 1 - cfg.varThresh;                  % 0.10 with default config
    d    = size(Xk, 2);
    mTh  = max(1, min(M - 1, d - 1));
    floor_ = max(4, 2*mTh + 2);                % fit's own statistical min
    segs = struct('fit', {}, 'size', {});      % field order fixed
    if size(Xk, 1) < max(4, floor_)
        return;                                % caller already gated, guard
    end
    stack = {Xk};
    while ~isempty(stack)
        X = stack{1};
        stack(1) = [];
        f = MRCoSEAManifoldFit(X, M, cfg);
        if f.valid && f.resid <= thr
            segs(end+1) = struct('fit', f, 'size', size(X, 1)); %#ok<AGROW>
        elseif f.m >= 1 && size(X, 1) >= 2*floor_ && ...
                numel(segs) + numel(stack) + 2 <= SEG_CAP
            z    = (X - f.mu) * f.V(:, 1);     % 1st principal coordinate
            [~, ord] = sort(z);                % stable: ties keep index
            half = floor(size(X, 1) / 2);
            stack{end+1} = X(ord(1:half), :);      %#ok<AGROW>
            stack{end+1} = X(ord(half+1:end), :);  %#ok<AGROW>
        else
            % not modelable and cannot be split further -> discard
        end
    end
end
