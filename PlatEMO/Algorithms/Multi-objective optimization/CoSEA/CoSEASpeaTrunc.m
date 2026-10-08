function keep = CoSEASpeaTrunc(Xn, nKeep)
%CoSEASpeaTrunc SPEA2-style iterative truncation in the decision space (v0.7).
%   keep = CoSEASpeaTrunc(Xn, nKeep) iteratively removes the individual whose
%   sorted (ascending) distance vector to the others is lexicographically
%   smallest, until nKeep individuals remain. Xn must be box-normalized.

    n    = size(Xn, 1);
    keep = true(1, n);
    while sum(keep) > nKeep
        idx = find(keep);
        m   = numel(idx);
        D   = pdist2(Xn(idx,:), Xn(idx,:));
        D(1:m+1:end) = inf;               % exclude self-distance
        Ds  = sort(D, 2);
        cand = 1:m;                       % lexicographic minimum
        for c = 1 : m
            mn   = min(Ds(cand, c));
            cand = cand(Ds(cand, c) == mn);
            if numel(cand) == 1
                break;
            end
        end
        keep(idx(cand(1))) = false;
    end
end
