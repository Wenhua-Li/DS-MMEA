function CrowdDis = CoSEACrowdingDec(Xn, FrontNo)
%CoSEACrowdingDec Axis-aligned crowding distance in the decision space (v0.6).
%   CrowdDis = CoSEACrowdingDec(Xn, FrontNo) replicates the old decision-space
%   crowding path kept for ablation: per-front axis-aligned crowding distance
%   with normalization by the population's own min/max per dimension (note:
%   this normalization drifts with generations, a known defect fixed in v0.7).

    n        = size(Xn, 1);
    CrowdDis = zeros(1, n);
    lo = min(Xn, [], 1);
    hi = max(Xn, [], 1);
    rngSpan = max(hi - lo, eps);
    for f = unique(FrontNo)
        idx = find(FrontNo == f);
        nF  = numel(idx);
        if nF <= 2
            CrowdDis(idx) = inf;
            continue;
        end
        Xf = (Xn(idx,:) - lo) ./ rngSpan;
        cd = zeros(nF, 1);
        for j = 1 : size(Xf, 2)
            [sv, so] = sort(Xf(:,j));
            cd(so(1))   = inf;
            cd(so(end)) = inf;
            for i = 2 : nF-1
                if isinf(cd(so(i)))
                    continue;
                end
                cd(so(i)) = cd(so(i)) + (sv(i+1) - sv(i-1));
            end
        end
        CrowdDis(idx) = cd.';
    end
end
