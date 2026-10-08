function [sel, elig] = MSCoSEAReport(Xn, Obj, Con, lab, N)
%MSCoSEAReport Report-layer mode-front filter + priority-tiered maximin.
%   [sel, elig] = MSCoSEAReport(Xn, Obj, Con, lab, N) assembles the index
%   set of the delivered population from the merged pool (design revision
%   2026-09-30, "prune at the report layer, never in the search"):
%
%     Xn  - [m x d] box-normalized decision vectors of the merged pool
%     Obj - [m x M] objectives of the merged pool
%     Con - [m x c] constraint values of the merged pool (or [])
%     lab - {1 x S} one row of cluster labels per reporting scale (the
%           wide family clustering of the merged pool, applied by the
%           caller; more scales generalize by union)
%     N   - size of the delivered set
%
%   Eligibility (who may enter the delivered set). The judgment is
%   MODE-RELATIVE wherever a mode carries information and POOL-RELATIVE
%   wherever it does not:
%     (a) a member of a REPORTABLE mode (a cluster with at least 3
%         members at some reporting scale) is eligible iff it is
%         nondominated WITHIN ITS OWN MODE -- the deliberately
%         mode-relative rule: on benchmarks whose truth includes local
%         Pareto fronts, globally dominated members are load-bearing
%         (71-problem deletion study, 2026-09-30);
%     (b) every other member -- unclustered (label 0) or in a trivial
%         cluster of 1 or 2 members -- is eligible iff it is on the
%         FRONT OF THE MERGED POOL. Within-mode sorting carries no
%         information there (two points are mutually nondominated
%         unless one dominates the other), so the pool-level front is
%         the only quality evidence available. This is the 2026-09-30
%         v4.2 tightening: the v4.1 rule let dominated noise members
%         and dominated pairs through unconditionally, which left
%         exactly the isolated poor solutions the filter exists to
%         remove (median 15-31 lone + 85-175 pair members on the
%         high-dimensional NMMF group, frozen-data diagnosis
%         diag_v41_residual). A real 2-member colony ON the pool front
%         is unaffected (it passes the pool test); only dominated
%         stragglers are pruned.
%
%   Fusion (who is actually delivered): greedy MAXIMIN in the normalized
%   decision space over the eligible pool. If fewer than N eligible
%   members exist, the same greedy rule CONTINUES over the remaining
%   pool in ascending pool-front order (the least-dominated leftovers
%   first, accumulated nearest-distance never reset), so the output is
%   always exactly N and the filter never silently disables itself.
%   Start point (rule unchanged from the original fusion): the eligible
%   member farthest from the overall centroid. Ties: smallest index
%   (deterministic). Zero parameters.
%
%   Returns sel (1 x N selected indices) and elig (1 x m eligibility).

    m    = size(Xn, 1);
    sel  = false(1, m);
    elig = false(1, m);

    %% Pool-level front numbers: quality evidence for noise/trivial
    %% members and the refill tier order
    if isempty(Con)
        gno = NDSort(Obj, m);
    else
        gno = NDSort(Obj, Con, m);
    end
    gno  = reshape(gno, 1, []);
    maxF = max(gno);

    %% Eligibility: mode-relative for reportable modes, pool-relative
    %% for noise and trivial clusters, union over reporting scales
    for s = 1 : numel(lab)
        ls  = reshape(lab{s}, 1, m);
        uC  = unique(ls(ls > 0));
        szs = zeros(1, max([uC, 0]));
        for c = uC(:).'
            szs(c) = nnz(ls == c);
        end
        eq   = (ls == 0);                    % noise: needs the pool front
        pos  = ls > 0;
        if any(pos)
            triv = false(1, m);
            triv(pos) = szs(ls(pos)) < 3;    % trivial clusters: same rule
            eq = eq | triv;
        end
        eq  = eq & (gno == 1);
        for c = uC(:).'
            if szs(c) < 3
                continue;                    % trivial: covered above
            end
            idx = find(ls == c);
            if isempty(Con)
                fno = NDSort(Obj(idx, :), numel(idx));
            else
                fno = NDSort(Obj(idx, :), Con(idx, :), numel(idx));
            end
            fno = reshape(fno, 1, []);
            eq(idx(fno == 1)) = true;
        end
        elig = elig | eq;
    end
    if ~any(elig)
        elig(:) = true;                      % degenerate guard: full pool
    end

    %% Priority-tiered greedy maximin (eligible first, then ascending
    %% pool fronts, accumulated dmin never reset)
    D = pdist2(Xn, Xn);
    D(1:m+1:end) = 0;
    d0 = reshape(vecnorm(Xn - mean(Xn, 1), 2, 2), 1, m);
    d0(~elig) = -inf;
    [~, j0] = max(d0);                  % start: eligible farthest (min idx)
    sel(j0) = true;
    dmin = D(j0, :);
    dmin(j0) = -inf;
    tier = 1;                            % refill walks the pool fronts
    for s = 2 : N
        cand = ~sel & elig;
        if ~any(cand)
            cand = ~sel & (gno == tier);
            while ~any(cand) && tier < maxF
                tier = tier + 1;
                cand = ~sel & (gno == tier);
            end
            if ~any(cand)
                cand = ~sel;            % safety net (fronts partition all)
            end
        end
        dj = dmin;
        dj(~cand) = -inf;
        [~, j] = max(dj);               % max picks the first max index
        sel(j) = true;
        dmin = min(dmin, D(j, :));
        dmin(sel) = -inf;
    end
    sel = find(sel);
end
