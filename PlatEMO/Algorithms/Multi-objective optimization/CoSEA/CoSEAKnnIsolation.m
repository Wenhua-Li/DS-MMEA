function Iso = CoSEAKnnIsolation(Xn, FrontNo)
%CoSEAKnnIsolation kNN isolation of each solution in the decision space (v0.7).
%   Iso = CoSEAKnnIsolation(Xn, FrontNo) returns the distance to the k-th
%   nearest neighbor within the same front, k = floor(sqrt(nF)), where nF is
%   the front size (SPEA2 convention). Single-member fronts receive Inf so
%   that newborn modes always win the tournament. Xn must be box-normalized.

    n   = size(Xn, 1);
    Iso = zeros(1, n);
    for f = unique(FrontNo)
        idx = find(FrontNo == f);
        nF  = numel(idx);
        if nF == 1
            Iso(idx) = inf;
            continue;
        end
        D = pdist2(Xn(idx,:), Xn(idx,:));
        D(1:nF+1:end) = inf;              % exclude self-distance
        k = min(max(1, floor(sqrt(nF))), nF-1);
        Ds = sort(D, 2);
        Iso(idx) = Ds(:,k).';
    end
end
